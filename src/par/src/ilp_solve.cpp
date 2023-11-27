#include <absl/types/span.h>
#include <ortools/sat/cp_model.h>
#include <ortools/util/sorted_interval_list.h>

#include <algorithm>
#include <csignal>
#include <cstddef>
#include <cstdint>
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
    void init(size_t num_partitions, CP::CpModelBuilder& cp_model)
    {
      m.resize(num_partitions);
      for (auto& row : m) {
        row.resize(num_partitions);
        for (auto& var : row) {
          var = cp_model.NewBoolVar();
        }
      }
    }
  };

  // Instantiate the solver
  CP::CpModelBuilder cp_model;

  auto tau = cp_model.NewIntVar(operations_research::Domain(0, in.max_delay));
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
      component_partition[component_idx][partition_idx] = cp_model.NewBoolVar();
    }
  }
  // Capacity constraints
  for (size_t part_idx = 0; part_idx < in.num_partitions; part_idx++) {
    auto res = CP::LinearExpr::WeightedSum(
        component_partition[part_idx],
        absl::Span<const int64_t>(component_utilization));
    cp_model.AddLessOrEqual(res, max_capacity[part_idx]);
  }

  // Instantiate edge variables for crossing constraints
  // We assume crossing only happens between the component and (replicated)
  // output buffer
  std::vector<Partition_Matrix_t> global_component_to_output(
      in.vertices.size());
  for (auto& edge : global_component_to_output) {
    edge.init(in.num_partitions, cp_model);
  }
  // Count the number of crossings
  auto crossig_domain = operations_research::Domain(0, in.vertices.size());
  std::vector<std::vector<CP::IntVar>> crossing_count(in.num_partitions);
  for (int src_part_idx = 0; src_part_idx < in.num_partitions; src_part_idx++) {
    crossing_count[src_part_idx].resize(in.num_partitions);
    for (int dst_part_idx = 0; dst_part_idx < src_part_idx; dst_part_idx++) {
      crossing_count[src_part_idx][dst_part_idx]
          = cp_model.NewIntVar(crossig_domain);
      operations_research::sat::LinearExpr crossing_count_expr;
      for (auto& edge : global_component_to_output) {
        crossing_count_expr += edge.m[src_part_idx][dst_part_idx];
        crossing_count_expr += edge.m[dst_part_idx][src_part_idx];
      }
      cp_model.AddLessOrEqual(
          crossing_count_expr,
          crossing_count[src_part_idx][dst_part_idx] * in.ifc_size);
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
      in.total_num_edges);
  {
    // Find the loops each edge belongs to
    std::vector<std::vector<int>> edge_loop_membership(in.total_num_edges);
    for (size_t loop_idx = 0; loop_idx < loops.size(); loop_idx++) {
      for (size_t vertex_iter_idx = 0; vertex_iter_idx < loops[loop_idx].size();
           vertex_iter_idx++) {
        edge_loop_membership[loops[loop_idx][vertex_iter_idx].vertex_idx]
            .push_back(loop_idx);
      }
    }
    // Populate provider variables
    for (size_t loop_idx = 0; loop_idx < in.total_num_edges; loop_idx++) {
      auto& provider_idxes = each_loop_edge_to_provider_idx[loop_idx];
      provider_idxes.resize(loops[loop_idx].size(), -1);
      for (const auto& edge : loops[loop_idx]) {
        auto& dst_vertex = in.vertices[edge.vertex_idx];
        // check whether any other loop needs this edge
        bool needed_by_other_loop = false;
        for (size_t fanin_idx = 0; fanin_idx < dst_vertex.fanins.size();
             fanin_idx++) {
          if (fanin_idx != edge.fanin_idx
              && !(
                  edge_loop_membership[dst_vertex.edge_id_start_idx + fanin_idx]
                      .empty())) {
            needed_by_other_loop = true;
            break;
          }
        }
        if (needed_by_other_loop) {
          auto cur_edge_idx = dst_vertex.edge_id_start_idx + edge.fanin_idx;
          auto& curr_edge_providers = edge_idx_to_provider_edges[cur_edge_idx];
          provider_idxes[cur_edge_idx] = curr_edge_providers.size();
          curr_edge_providers.emplace_back(in.num_partitions, cp_model);
        }
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
    auto output_port_var = make_var_vec();
    auto component_var = make_var_vec();
    std::vector<CP::IntVar> within_loop_delay(loop_size);
    std::generate(within_loop_delay.begin(),
                  within_loop_delay.end(),
                  [&cp_model, delay_domain]() {
                    return cp_model.NewIntVar(delay_domain);
                  });

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
      auto& cur_component_var = component_var[vertex_idx];
      // local imply global
      // component partition
      for (int part_idx = 0; part_idx < num_partitions; part_idx++) {
        cp_model.AddImplication(
            cur_component_var[part_idx],
            component_partition[dst_vertex.foreign_id.cluster_id][part_idx]);
      }
      // component at partition a and output at partition b
      for (int part_a_idx = 0; part_a_idx < num_partitions; part_a_idx++) {
        for (int part_b_idx = 0; part_b_idx < num_partitions; part_b_idx++) {
          CP::BoolVar local_edge[] = {cur_component_var[part_a_idx],
                                      output_port_var[vertex_idx][part_b_idx]};
          CP::BoolVar global_edge[]
              = {global_component_to_output[dst_vertex.foreign_id.cluster_id]
                     .m[part_a_idx][part_b_idx]};
          cp_model.AddImplication(local_edge, global_edge);

          // Within loop delay
          cp_model
              .AddLessOrEqual(
                  in.fixed_crossing_delay
                      + in.proportional_crossing_delay
                            * crossing_count[part_a_idx][part_b_idx],
                  within_loop_delay[vertex_idx])
              .OnlyEnforceIf(local_edge);
        }
      }

      // A component in partition a implies all its fanin terminates at
      // partition a
      for (int fanin_idx = 0; fanin_idx < dst_vertex.fanins.size();
           fanin_idx++) {
        if (fanin_idx == edge.fanin_idx) {
          // Loop within current loop
          for (int part_idx = 0; part_idx < num_partitions; part_idx++) {
            cp_model.AddImplication(cur_component_var[part_idx],
                                    output_port_var[vertex_idx - 1][part_idx]);
          }
          continue;
        }
        auto cur_edge_idx = dst_vertex.edge_id_start_idx + fanin_idx;
        auto& curr_edge_providers = edge_idx_to_provider_edges[cur_edge_idx];
        std::vector<CP::LinearExpr> expr(num_partitions);
        if (curr_edge_providers.empty()) {
          raise(SIGTRAP);
        }
        for (auto& provider : curr_edge_providers) {
          for (auto& edges : provider.m) {
            for (int part_idx = 0; part_idx < num_partitions; part_idx++) {
              expr[part_idx] += edges[part_idx];
            }
          }
        }
        for (int part_idx = 0; part_idx < num_partitions; part_idx++) {
          cp_model.AddLessOrEqual(cur_component_var[part_idx], expr[part_idx]);
        }
      }

      // Within Loop Delay
      cp_model.AddLessOrEqual(
          CP::LinearExpr::Sum(within_loop_delay) + loop_comb_delay,
          tau * loop_reg_cnt);

      // Enforce source presence for edges feeding into other loops
      //These edges are considered as running in parallel to the component to output buffer edges
      if (loop_edge_to_provider_idx[vertex_idx] != -1) {
        auto cur_edge_idx = dst_vertex.edge_id_start_idx + edge.fanin_idx;
        auto& to_other_loop_edge
            = edge_idx_to_provider_edges[cur_edge_idx]
                                        [loop_edge_to_provider_idx[vertex_idx]];
        auto cross_loop_edge_delay = cp_model.NewIntVar(delay_domain);

        for (int src_part_idx = 0; src_part_idx < num_partitions;
             src_part_idx++) {
          for (int dst_part_idx = 0; dst_part_idx < num_partitions;
               dst_part_idx++) {
            auto& dst_var = to_other_loop_edge.m[src_part_idx][dst_part_idx];
            // The edge (dst_var) start from src_part_idx means the driver component is at src_part_idx
            cp_model.AddImplication(
                dst_var, component_var[vertex_idx - 1][src_part_idx]);
            //Add the corresponding edge crossing count
            auto driver_idx=dst_vertex.fanins[edge.fanin_idx].driver_id;
            cp_model.AddImplication(
                dst_var, global_component_to_output[driver_idx].m[src_part_idx][dst_part_idx]);
            cp_model
                .AddLessOrEqual(crossing_count[src_part_idx][dst_part_idx]
                                        * in.proportional_crossing_delay
                                    + in.fixed_crossing_delay,
                                cross_loop_edge_delay)
                .OnlyEnforceIf(dst_var);
          }
        }
        //Cross loop delay
        //Include all edges except the one parallel to the current edge reaching other loops
        auto prev_vertex_idx = vertex_idx - 1;
        if (prev_vertex_idx ==-1) {
          prev_vertex_idx = loop_size-1;
        }
        CP::LinearExpr cross_loop_delay_expr=loop_comb_delay+cross_loop_edge_delay;
        for (int vertex_idx = 0; vertex_idx < loop_size; vertex_idx++) {
          if (vertex_idx == prev_vertex_idx) {
            continue;
          }
          cross_loop_delay_expr += within_loop_delay[vertex_idx];
        }
        cp_model.AddLessOrEqual(cross_loop_delay_expr, tau * loop_reg_cnt);
      }
    }
  }

  cp_model.Minimize(tau);
  
}
