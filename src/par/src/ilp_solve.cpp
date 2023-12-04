#include <spdlog/fmt/bundled/core.h>
#include <spdlog/fmt/bundled/format.h>
#include <spdlog/fmt/bundled/color.h>
#include <spdlog/fmt/bundled/ostream.h>
#include <absl/types/span.h>
#include <ortools/sat/cp_model.h>
#include <ortools/sat/cp_model.pb.h>
#include <ortools/sat/cp_model_solver.h>
#include <ortools/sat/sat_parameters.pb.h>
#include <ortools/util/sorted_interval_list.h>


#include <algorithm>
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "replicate_retime_partition.hpp"
namespace CP = operations_research::sat;
template <typename T>
struct Circular_Vector : public std::vector<T>
{
  using std::vector<T>::vector;
  T& operator[](size_t idx)
  {
    if (idx == -1) {
      return this->back();
    }
    if (idx >= this->size()) {
      raise(SIGTRAP);
    }
    return std::vector<T>::operator[](idx);
  }
};
void par::RRP::ilp_part(const Graph& in,
                        Partition_Solution_t& sol_out,
                        const std::vector<std::vector<Edge>>& loops,
                        const std::vector<size_t>& max_capacity,
                        const std::vector<int64_t>& component_utilization)
{
  struct Partition_Matrix_t
  {
    std::vector<std::vector<CP::BoolVar>> m;
    Partition_Matrix_t() = default;
    Partition_Matrix_t(int num_partitions, CP::CpModelBuilder& cp_model,const std::string& name = ""){
      init(num_partitions, cp_model, name);
    }
    void init(int num_partitions, CP::CpModelBuilder& cp_model,const std::string& name = "")
    {
      m.resize(num_partitions);
      for (int src_idx=0;src_idx<num_partitions;src_idx++) {
        auto &row = m[src_idx];
        row.resize(num_partitions);
        for (int dst_idx=0;dst_idx<num_partitions;dst_idx++) {
          auto &var = row[dst_idx];
          var = cp_model.NewBoolVar();
          if(!name.empty()){
            var.WithName(fmt::format("{}_p{}_p{}", name, src_idx, dst_idx));
          }
        }
      }
    }
  };

  // Instantiate the solver
  CP::CpModelBuilder cp_model;

  auto tau = cp_model.NewIntVar(operations_research::Domain(0, in.max_delay)).WithName("tau");
  // Instatiate global components variables for capacity constraints
  // The component has size, the output is a size-less buffer to allow a batch
  // of fanout going to loads on one partition to have only one crossing via its
  // replication
  std::vector<std::vector<CP::BoolVar>> component_partition(
      in.max_component_cnt);
  for (size_t component_idx = 0; component_idx < in.max_component_cnt;
       component_idx++) {
    component_partition[component_idx].resize(in.num_partitions);
    for (size_t partition_idx = 0; partition_idx < in.num_partitions;
         partition_idx++) {
      component_partition[component_idx][partition_idx] = cp_model.NewBoolVar().WithName(
          fmt::format("c{}_p_{}", component_idx, partition_idx));
    }
  }
  // Capacity constraints
  for (size_t part_idx = 0; part_idx < in.num_partitions; part_idx++) {
    CP::LinearExpr expr;
    for(size_t component_idx = 0; component_idx < in.max_component_cnt;
        component_idx++){
      expr += component_partition[component_idx][part_idx]*component_utilization[component_idx];
    }
    
    cp_model.AddLessOrEqual(expr, max_capacity[part_idx]).WithName(
        fmt::format("capacity_p{}", part_idx));
  }

  // Instantiate edge variables for crossing constraints
  // We assume crossing only happens between the component and (replicated)
  // output buffer
  std::vector<Partition_Matrix_t> global_component_to_output(
      in.vertices.size());
  for (size_t vertex_idx = 0; vertex_idx < in.vertices.size(); vertex_idx++) {
    auto& vertex_out = global_component_to_output[vertex_idx];
    vertex_out.init(in.num_partitions, cp_model,fmt::format("vout_{}_global_crossing",vertex_idx));
    //Each output partition only needs to be sourced from one partition
    for (int dst_part_idx = 0; dst_part_idx < in.num_partitions; dst_part_idx++) {
      CP::LinearExpr expr;
      for (int src_part_idx = 0; src_part_idx < in.num_partitions; src_part_idx++) {
        expr += vertex_out.m[src_part_idx][dst_part_idx];
      }
      cp_model.AddLessOrEqual(expr, 1).WithName(
          fmt::format("vout_{}_p{}_src", vertex_idx, dst_part_idx));
    }
  }
  // Count the number of crossings
  auto crossig_domain = operations_research::Domain(0, in.vertices.size());
  std::vector<std::vector<CP::IntVar>> crossing_count(in.num_partitions);
  for (int src_part_idx = 0; src_part_idx < in.num_partitions; src_part_idx++) {
    crossing_count[src_part_idx].resize(in.num_partitions);
    for (int dst_part_idx = 0; dst_part_idx < src_part_idx; dst_part_idx++) {
      crossing_count[src_part_idx][dst_part_idx]
          = cp_model.NewIntVar(crossig_domain).WithName(fmt::format("crossing_p{}_p{}", src_part_idx, dst_part_idx));
      operations_research::sat::LinearExpr crossing_count_expr;
      for (auto& edge : global_component_to_output) {
        crossing_count_expr += edge.m[src_part_idx][dst_part_idx];
        crossing_count_expr += edge.m[dst_part_idx][src_part_idx];
      }
      cp_model.AddLessOrEqual(
          crossing_count_expr,
          crossing_count[src_part_idx][dst_part_idx] * in.ifc_size).WithName(fmt::format("crossing_p{}_p{}_sum", src_part_idx, dst_part_idx));;
    }
  }
  for (int src_part_idx = 0; src_part_idx < in.num_partitions; src_part_idx++) {
    for (int dst_part_idx = src_part_idx + 1; dst_part_idx < in.num_partitions;
         dst_part_idx++) {
      crossing_count[src_part_idx][dst_part_idx]
          = crossing_count[dst_part_idx][src_part_idx];
    }
  }

  // Declare Variables for edges between loops
  std::vector<std::vector<Partition_Matrix_t>> edge_idx_to_provider_edges(
      in.total_num_edges);
  std::vector<std::vector<int>> each_loop_edge_to_provider_idx(
      loops.size());
  {
    // Find the loops each edge belongs to
    /*std::vector<std::vector<int>> edge_loop_membership(in.total_num_edges);
    for (size_t loop_idx = 0; loop_idx < loops.size(); loop_idx++) {
      for (size_t vertex_iter_idx = 0; vertex_iter_idx < loops[loop_idx].size();
           vertex_iter_idx++) {
        edge_loop_membership[loops[loop_idx][vertex_iter_idx].vertex_idx]
            .push_back(loop_idx);
      }
    }*/
    // Populate provider variables
    for (size_t loop_idx = 0; loop_idx < loops.size(); loop_idx++) {
      const auto& loop = loops[loop_idx];
      printf( "Loop %zu:", loop_idx);
      for (const auto& vertex : loop) {
        printf( " %zu,", vertex.vertex_idx);
      }
      printf( "\n");
      auto& provider_idxes = each_loop_edge_to_provider_idx[loop_idx];
      provider_idxes.resize(loops[loop_idx].size(), -1);
      int edge_idx=0;
      for (const auto& edge : loops[loop_idx]) {
        auto& dst_vertex = in.vertices[edge.vertex_idx];
        // check whether any other loop needs this edge
        /*bool needed_by_other_loop = false;
        for (size_t fanin_idx = 0; fanin_idx < dst_vertex.fanins.size();
             fanin_idx++) {
          if (fanin_idx != edge.fanin_idx) {
            needed_by_other_loop = true;
            break;
          }
        }
        if (needed_by_other_loop) {*/
          auto cur_edge_idx = dst_vertex.edge_id_start_idx + edge.fanin_idx;
          auto& curr_edge_providers = edge_idx_to_provider_edges[cur_edge_idx];
          if(curr_edge_providers.size()>loop_idx){
            raise(SIGTRAP);
          }
          provider_idxes[edge_idx] = curr_edge_providers.size();
          curr_edge_providers.emplace_back(in.num_partitions, cp_model, fmt::format("l{}_v{}tov{}_provider", 
            loop_idx, dst_vertex.fanins[edge.fanin_idx].driver_id, edge.vertex_idx));
        //}
        edge_idx++;
      }
    }
  }

  operations_research::Domain delay_domain(0, in.max_crossing);
  // Variables for each loop:
  for (size_t loop_idx = 0; loop_idx < loops.size(); loop_idx++) {
    const auto& loop = loops[loop_idx];
    const auto& loop_edge_to_provider_idx
        = each_loop_edge_to_provider_idx[loop_idx];
    auto num_partitions = in.num_partitions;
    auto loop_size = loop.size();
    auto make_var_vec = [&cp_model, num_partitions, loop_size]() {
      std::vector<std::vector<CP::BoolVar>> res(loop_size);
      std::generate(res.begin(), res.end(), [&cp_model, num_partitions]() {
        std::vector<CP::BoolVar> res(num_partitions);
        std::generate(res.begin(), res.end(), [&cp_model]() {
          return cp_model.NewBoolVar();
        });
        return res;
      });
      return res;
    };
    std::vector<std::vector<CP::BoolVar>> output_port_var(loop_size);
    std::vector<CP::IntVar> within_loop_delay(loop_size);
    for(int vertex_idx=0; vertex_idx<loop_size; vertex_idx++) {
      output_port_var[vertex_idx].resize(in.num_partitions);
      for(int part_idx=0; part_idx<in.num_partitions; part_idx++) {
        auto & vert_var = output_port_var[vertex_idx][part_idx];
        vert_var = cp_model.NewBoolVar();
        vert_var.WithName(fmt::format("vertex_l{}_v{}_p{}", loop_idx, loop[vertex_idx].vertex_idx, part_idx));
        auto & edge_delay_var=within_loop_delay[vertex_idx];
        edge_delay_var = cp_model.NewIntVar(delay_domain);
        auto next_vertex_idx = (vertex_idx+1)%loop_size;
        edge_delay_var.WithName(fmt::format("within_loop_l{}_v{}tov{}_delay", loop_idx, loop[vertex_idx].vertex_idx, loop[next_vertex_idx].vertex_idx));

      }
    }
    // Orignal delay of the loop without cuts
    auto loop_comb_delay = 0;
    auto loop_reg_cnt = 0;
    for (const auto& edge : loop) {
      auto& dst_vertex = in.vertices[edge.vertex_idx];
      auto& fanin_edge = dst_vertex.fanins[edge.fanin_idx];
      loop_reg_cnt += fanin_edge.reg_cnt;
      loop_comb_delay += fanin_edge.comb_delay;
    }

    for (int vertex_idx = 0; vertex_idx < loop_size; vertex_idx++) {
      auto& edge = loop[vertex_idx];
      auto& dst_vertex = in.vertices[edge.vertex_idx];
      auto& cur_output_port_var = output_port_var[vertex_idx];
      auto next_idx=(vertex_idx + 1) % loop_size;
      auto& next_output_port_var= output_port_var[next_idx];
      //auto& cur_component_var = component_var[vertex_idx];
      // local imply global
      // component partition
      for (int part_idx = 0; part_idx < num_partitions; part_idx++) {
        const auto& res=cp_model.AddImplication(
            cur_output_port_var[part_idx],
            component_partition[dst_vertex.foreign_id.cluster_id][part_idx]).WithName(fmt::format("v{}_l{}_p{}_util_to_g", dst_vertex.id, loop_idx, part_idx));
        if(part_idx==0){
          fprintf(stderr, "%s:%s=>%s\n",res.Name().c_str(),cur_output_port_var[part_idx].Name().c_str(),component_partition[dst_vertex.foreign_id.cluster_id][part_idx].Name().c_str());
        }
      }
      cp_model.AddAtLeastOne(cur_output_port_var).WithName(fmt::format("v{}_l{}_present", dst_vertex.id, loop_idx));
      // cur_output at partition a and next_output at partition b
      for (int part_a_idx = 0; part_a_idx < num_partitions; part_a_idx++) {
        for (int part_b_idx = 0; part_b_idx < num_partitions; part_b_idx++) {
          CP::BoolVar local_edge[] = {cur_output_port_var[part_a_idx],
                                      next_output_port_var[part_b_idx]};
          CP::BoolVar global_edge[]
              = {global_component_to_output[dst_vertex.foreign_id.cluster_id]
                     .m[part_a_idx][part_b_idx]};
          const auto& cons=cp_model.AddImplication(local_edge, global_edge).WithName(fmt::format("within_loop_{}_v{}tov{}_p{}p{}_global_edgeq",loop_idx, dst_vertex.id,loop[next_idx].vertex_idx, part_a_idx, part_b_idx));
          if(part_a_idx==0&&part_b_idx==0){
            fprintf(stderr, "%s: %s&&%s=>%s\n",cons.Name().c_str(),local_edge[0].Name().c_str(),local_edge[1].Name().c_str(),global_edge[0].Name().c_str());
          }
          if(part_a_idx!=part_b_idx){
          // Within loop delay
          cp_model
              .AddLessOrEqual(
                  in.fixed_crossing_delay
                      + in.proportional_crossing_delay
                            * crossing_count[part_a_idx][part_b_idx],
                  within_loop_delay[vertex_idx])
              .OnlyEnforceIf(local_edge).WithName(fmt::format("within_loop_{}_v{}tov{}_p{}p{}_delay_cstr",loop_idx, dst_vertex.id,loop[next_idx].vertex_idx, part_a_idx, part_b_idx));
          }
        }
      }

      // A component in partition a implies all its fanin terminates at
      // partition a
      for (int fanin_idx = 0; fanin_idx < dst_vertex.fanins.size();
           fanin_idx++) {
        if (fanin_idx == edge.fanin_idx) {
          // crossing on edge within the current loop is already accounted for in the previous step, so no additional constraint for within loop crossing delay 
          continue;
        }
        //Need fanins from other loops
        auto cur_edge_idx = dst_vertex.edge_id_start_idx + fanin_idx;
        auto& curr_edge_providers = edge_idx_to_provider_edges[cur_edge_idx];
        std::vector<CP::LinearExpr> expr(num_partitions);
        if (curr_edge_providers.empty()) {
          fmt::print(stderr,fmt::emphasis::bold | fg(fmt::color::red),
            "in loop {}, no provider for edge {} to {}\n",
            loop_idx,dst_vertex.fanins[fanin_idx].driver_id,dst_vertex.id);
          continue;
        }
        std::string provider_str;
        //Only need one of the providers to be true
        for (auto& provider : curr_edge_providers) {
          //And provide from one of the source partitions
          for (auto& edges : provider.m) {
            for (int part_idx = 0; part_idx < num_partitions; part_idx++) {
              expr[part_idx] += edges[part_idx];
              if (part_idx == 0) {
                provider_str += edges[part_idx].Name() + "+";
              }              
            }
          }
        }
        auto& cur_fanin=dst_vertex.fanins[fanin_idx];
        for (int part_idx = 0; part_idx < num_partitions; part_idx++) {
          const auto& cons=cp_model.AddLessOrEqual(cur_output_port_var[part_idx], expr[part_idx]).WithName(
            fmt::format("l{}_v{}tov{}_p{}_other_loops", loop_idx,cur_fanin.driver_id,dst_vertex.id, part_idx));
          if(part_idx==0){
            fprintf(stderr, "%s: %s<=%s\n",cons.Name().c_str(),cur_output_port_var[part_idx].Name().c_str(),provider_str.c_str());
          }
        }
      }

      // Within Loop Delay
      cp_model.AddLessOrEqual(
          CP::LinearExpr::Sum(within_loop_delay) + loop_comb_delay,
          tau * loop_reg_cnt).WithName(fmt::format("l{}_inside_delay", loop_idx));

      // Enforce source presence for edges feeding into other loops
      //These edges are considered as running in parallel to the previous port to current port output buffer edges
      if (loop_edge_to_provider_idx[vertex_idx] != -1) {
        auto cur_edge_idx = dst_vertex.edge_id_start_idx + edge.fanin_idx;
        auto& to_other_loop_edge
            = edge_idx_to_provider_edges[cur_edge_idx]
                                        [loop_edge_to_provider_idx[vertex_idx]];
        auto cross_loop_edge_delay = cp_model.NewIntVar(delay_domain).WithName(
            fmt::format("l{}_v{}tov{}_cross_loop_delay", loop_idx,dst_vertex.fanins[edge.fanin_idx].driver_id ,dst_vertex.id));
        
        auto prev_vertex_idx = vertex_idx - 1;
        if (prev_vertex_idx ==-1) {
          prev_vertex_idx = loop_size-1;
        }
        
        for (int src_part_idx = 0; src_part_idx < num_partitions;
             src_part_idx++) {
          for (int dst_part_idx = 0; dst_part_idx < num_partitions;
               dst_part_idx++) {
            auto& dst_var = to_other_loop_edge.m[src_part_idx][dst_part_idx];
            // The edge (dst_var) start from src_part_idx means the driver component is at src_part_idx
            const auto& cons1=cp_model.AddImplication(
                dst_var, output_port_var[prev_vertex_idx][src_part_idx]).WithName(
                fmt::format("l{}_v{}tov{}_p{}p{}_require_src_present", loop_idx,dst_vertex.fanins[edge.fanin_idx].driver_id ,dst_vertex.id, src_part_idx, dst_part_idx)
                );
            if(src_part_idx==0&&dst_part_idx==0){
              fprintf(stderr, "%s: %s=>%s\n",cons1.Name().c_str(),dst_var.Name().c_str(),output_port_var[prev_vertex_idx][src_part_idx].Name().c_str());
            }
            //Add the corresponding edge crossing count
            auto driver_idx=dst_vertex.fanins[edge.fanin_idx].driver_id;
            const auto& cons2=cp_model.AddImplication(
                dst_var, global_component_to_output[driver_idx].m[src_part_idx][dst_part_idx]).WithName(
                fmt::format("l{}_v{}tov{}_p{}p{}_require_global_crossing", loop_idx,dst_vertex.fanins[edge.fanin_idx].driver_id ,dst_vertex.id, src_part_idx, dst_part_idx)
                );
            if(src_part_idx==0&&dst_part_idx==0){
              fprintf(stderr, "%s: %s=>%s\n",cons2.Name().c_str(),dst_var.Name().c_str(),global_component_to_output[driver_idx].m[src_part_idx][dst_part_idx].Name().c_str());
            }
            if(src_part_idx!=dst_part_idx){
            cp_model
                .AddLessOrEqual(crossing_count[src_part_idx][dst_part_idx]
                                        * in.proportional_crossing_delay
                                    + in.fixed_crossing_delay,
                                cross_loop_edge_delay)
                .OnlyEnforceIf(dst_var).WithName(
                fmt::format("l{}_v{}tov{}_p{}p{}_cross_loop_delay", loop_idx,dst_vertex.fanins[edge.fanin_idx].driver_id ,dst_vertex.id, src_part_idx, dst_part_idx)
                );
            }
          }
        }
        //Cross loop delay
        //Include all edges except the one parallel to the current edge reaching other loops
        std::stringstream cross_loop_delay_ss;
        cross_loop_delay_ss<<loop_comb_delay<<'+';
        CP::LinearExpr cross_loop_delay_expr=loop_comb_delay+cross_loop_edge_delay;
        for (int vertex_idx = 0; vertex_idx < loop_size; vertex_idx++) {
          if (vertex_idx == prev_vertex_idx) {
            continue;
          }
          cross_loop_delay_expr += within_loop_delay[vertex_idx];
          cross_loop_delay_ss<<within_loop_delay[vertex_idx].Name()<<'+';
        }
        const auto& cons=cp_model.AddLessOrEqual(cross_loop_delay_expr, tau * loop_reg_cnt).WithName(
            fmt::format("l{}_v{}tov{}_cross_loop_delay", loop_idx,dst_vertex.fanins[edge.fanin_idx].driver_id ,dst_vertex.id)
            );
        fprintf(stderr, "%s:%s<=tau*%d\n",cons.Name().c_str(),cross_loop_delay_ss.str().c_str(),loop_reg_cnt);
      }
    }
  }

  cp_model.Minimize(tau);

  CP::SatParameters params;
  params.set_max_time_in_seconds(600);
  //params.set_log_search_progress(true);
  auto built_model = cp_model.Build();
  auto fh=fopen("model_constraints", "w");
  fprintf(fh, "%s\n", built_model.DebugString().c_str());
  fclose(fh);
  auto response = CP::SolveWithParameters(built_model,params);
  sol_out.component_to_output_solutions.clear();
  sol_out.component_to_output_solutions.resize(in.vertices.size());
  if (response.status() == CP::OPTIMAL || response.status() == CP::FEASIBLE) {
    sol_out.clock_period = response.objective_value();
    for (size_t edge_idx = 0; edge_idx < global_component_to_output.size();
         edge_idx++) {
      auto& edge_var = global_component_to_output[edge_idx];
      for (int src_part_idx = 0; src_part_idx < in.num_partitions;
           src_part_idx++) {
        for (int dst_part_idx = 0; dst_part_idx < in.num_partitions;
             dst_part_idx++) {
          if (CP::SolutionBooleanValue(
                  response, edge_var.m[src_part_idx][dst_part_idx])) {
            sol_out.component_to_output_solutions[edge_idx]
                .partitions.emplace_back(src_part_idx, dst_part_idx);
          }
        }
      }
    }
    sol_out.crossing_count.resize(in.num_partitions);
    for (int src_part_idx = 0; src_part_idx < in.num_partitions;
         src_part_idx++) {
      sol_out.crossing_count[src_part_idx].resize(in.num_partitions);
      for (int dst_part_idx = 0; dst_part_idx < in.num_partitions;
           dst_part_idx++) {
            sol_out.crossing_count[src_part_idx][dst_part_idx]=CP::SolutionIntegerValue(response,crossing_count[src_part_idx][dst_part_idx]);
      }
    }
  auto fh=fopen("solve_all_var", "w");
  for(int idx=0;idx<response.solution_size();idx++){
    if(response.solution(idx)){
      fprintf(fh, "%s=%ld\n",built_model.variables(idx).name().c_str(),response.solution(idx));
    }
  }
  fclose(fh);
  } else {
    std::cout << response.solve_log();
    std::cout << "No optimal solution found" << std::endl;
  }
}
