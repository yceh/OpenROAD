#include <spdlog/fmt/bundled/ostream.h>
#include <spdlog/fmt/bundled/color.h>
#include <algorithm>
#include <boost/container_hash/extensions.hpp>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <unordered_set>
#include <vector>

#include "replicate_retime_partition.hpp"
#include "signal.h"
// Given a partition solution, construct a graph after replication and account
// for crossing delays. Then, find loops whose delay> # of cycles* clock period.
namespace Find_Violating_Loops {
struct Vertex_Attributes
{
  size_t ori_vertex_idx;  // Multiple vertices in the replicated graph may
                          // correspond to the same vertex in the original
                          // graph, they will have different partition ids.
  int partition_id;
  int time_offset;
  Vertex_Attributes(size_t ori_vertex_idx, int partition_id)
      : ori_vertex_idx(ori_vertex_idx),
        partition_id(partition_id),
        time_offset(0){};
};

struct Retiming_Graph
{
  std::vector<Vertex_Attributes> attributes;
  par::RRP::Graph g;
};
void construct_graph_after_replication(
    par::RRP::Partition_Solution_t& partition_solution,
    const par::RRP::Graph& in,
    Retiming_Graph& g_out)
{
  class Idx_Map_t
  {
    std::vector<size_t> new_idxes;
    int num_partitions;
    Retiming_Graph& g_out;

   public:
    Idx_Map_t(const par::RRP::Graph& in, Retiming_Graph& g_out)
        : new_idxes(in.vertices.size() * in.num_partitions, SIZE_MAX),
          num_partitions(in.num_partitions),
          g_out(g_out)
    {
    }
    size_t operator()(size_t ori_idx, int partition_id) const
    {
      return new_idxes[ori_idx * num_partitions + partition_id];
    }
    bool add_vertex(size_t ori_idx, int partition_id)
    {
      auto idx = ori_idx * num_partitions + partition_id;
      if (new_idxes[idx] == SIZE_MAX) {
        new_idxes[idx] = g_out.g.vertices.size();
        g_out.g.vertices.emplace_back();
        g_out.attributes.emplace_back(ori_idx, partition_id);
        g_out.g.vertices.back().id = new_idxes[idx];
        return true;
      }
      return false;
    }
  };
  Idx_Map_t idx_map(in, g_out);
  auto& out_vertices = g_out.g.vertices;
  out_vertices.reserve(in.num_partitions * in.vertices.size());
  g_out.attributes.reserve(in.num_partitions * in.vertices.size());
  for (size_t out_idx = 0; out_idx < in.vertices.size(); ++out_idx) {
    for (const auto& partition_crossings :
         partition_solution.component_to_output_solutions[out_idx].partitions) {
      idx_map.add_vertex(out_idx, partition_crossings.second);
    }
  }
  //There may be edges not includes in the loops in the last iteration, need to route them
  for (size_t out_idx = 0; out_idx < in.vertices.size(); ++out_idx) {
    for (const auto& partition_crossings :
         partition_solution.component_to_output_solutions[out_idx].partitions) {
      const auto& ori_fanin = in.vertices[out_idx].fanins;
      for (const auto& fanin: ori_fanin) {
        if(idx_map.add_vertex(fanin.driver_id, partition_crossings.first)){
          fmt::print(stderr,fmt::emphasis::bold | fg(fmt::color::red),
            "driver {} of {} in part {} not in {}\n",
            fanin.driver_id, out_idx, partition_crossings.second,partition_crossings.first);
          auto& fanin_parts=partition_solution.component_to_output_solutions[fanin.driver_id].partitions;
          //just pick the first one, edge may be added in the next iteration
          auto src_part=fanin_parts.begin()->first;
          fanin_parts.emplace_back(src_part,partition_crossings.first);
          partition_solution.crossing_count[src_part][partition_crossings.first]++;
        }
      }
      }
  }

  for (size_t out_idx = 0; out_idx < in.vertices.size(); ++out_idx) {
    for (const auto& partition_crossings :
         partition_solution.component_to_output_solutions[out_idx].partitions) {
      auto dst_idx = idx_map(out_idx, partition_crossings.second);
      auto& dst_vertice = out_vertices[dst_idx];
      if (!dst_vertice.fanins.empty()) {
        raise(SIGTRAP);
      }
      auto crossing_delay
          = partition_crossings.first == partition_crossings.second
                ? 0
                : (in.fixed_crossing_delay
                   + in.proportional_crossing_delay
                         * partition_solution
                               .crossing_count[partition_crossings.first]
                                              [partition_crossings.second]);

      const auto& ori_fanin = in.vertices[out_idx].fanins;
      dst_vertice.fanins.resize(ori_fanin.size());
      for (size_t fanin_idx = 0; fanin_idx < ori_fanin.size(); fanin_idx++) {
        auto& cur_fanin_out = dst_vertice.fanins[fanin_idx];
        auto& cur_fanin_ori = ori_fanin[fanin_idx];

        // If the comb delay containt multiple clock cycles, subtract them from
        // the register count (it has to if timing is satisfied)
        auto new_comb_delay = cur_fanin_ori.comb_delay + crossing_delay;
        auto reg_cnt_substracted
            = std::min(cur_fanin_ori.reg_cnt,
                       new_comb_delay / partition_solution.clock_period);
        new_comb_delay -= reg_cnt_substracted * partition_solution.clock_period;
        cur_fanin_out.comb_delay = new_comb_delay;
        cur_fanin_out.reg_cnt = cur_fanin_ori.reg_cnt - reg_cnt_substracted;
        cur_fanin_out.driver_id
            = idx_map(cur_fanin_ori.driver_id, partition_crossings.first);
        if(cur_fanin_out.driver_id == SIZE_MAX){
          raise(SIGTRAP);
        }
      }
    }
  }
  par::RRP::populate_fanout(g_out.g);
}

}  // namespace Find_Violating_Loops

bool par::RRP::find_violating_loops(
    par::RRP::Partition_Solution_t& partition_solution,
    const par::RRP::Graph& in,
    std::vector<std::vector<par::RRP::Edge>>& violating_loops)
{
  Find_Violating_Loops::Retiming_Graph g_out;
  Find_Violating_Loops::construct_graph_after_replication(
      partition_solution, in, g_out);
  std::vector<std::vector<par::RRP::Edge>> loops;
  par::RRP::find_negtive_slack_loops(
      g_out.g, partition_solution.clock_period, loops);

  // Translate the loops back to the original graph
  typedef std::vector<size_t> Loop_t;
  struct Loop_Container
  {
    std::unordered_set<Loop_t, boost::hash<Loop_t>> loops;
    bool added_any = false;
    void operator()(Loop_t::iterator begin, Loop_t::iterator end)
    {
      Loop_t loop(begin, end);
      if (loop.size() < 2) {
        raise(SIGTRAP);
      }
      auto min_element_iter = std::min_element(loop.begin(), loop.end());
      std::rotate(loop.begin(), min_element_iter, loop.end());
      added_any |= loops.emplace(std::move(loop)).second;
    }
  };
  Loop_Container loop_container;
  loop_container.added_any = false;
  for (const auto& loop : loops) {
    // Decompose loop
    // Translate to original vertex id
    Loop_t ori_vertex_idx(loop.size());
    fprintf(stderr, "===Violating loop==========================\n");
    int total_reg=0,total_delay=0;
    for (size_t i = 0; i < loop.size(); i++) {
      auto& attr=g_out.attributes[loop[i].vertex_idx];
      ori_vertex_idx[i] = attr.ori_vertex_idx;
      auto& fanin=g_out.g.vertices[loop[i].vertex_idx].fanins[loop[i].fanin_idx];
      total_delay+=fanin.comb_delay;
      total_reg+=fanin.reg_cnt;
      fprintf(stderr, "\t vertex %zu in part %d,delay %d,reg %d\n",
              attr.ori_vertex_idx, attr.partition_id,fanin.comb_delay,fanin.reg_cnt);
    }
    fprintf(stderr, "total delay %d, total reg %d, slack %d\n\n",total_delay,total_reg,partition_solution.clock_period*total_reg-total_delay);
    // Excise sub loops
    bool found = false;
    do {
      found = false;
      std::unordered_map<size_t, Loop_t::iterator> prev_vertex_idx(
          ori_vertex_idx.size());
      for (auto iter = ori_vertex_idx.begin(); iter != ori_vertex_idx.end();
           iter++) {
        auto res = prev_vertex_idx.emplace(*iter, iter);
        if (!res.second) {
          // Found loop
          found = true;
          auto start_iter = res.first->second;
          loop_container(start_iter, iter);
          ori_vertex_idx.erase(start_iter, iter);
          break;
        }
      }
    } while (found);
    if (!ori_vertex_idx.empty()) {
      loop_container(ori_vertex_idx.begin(), ori_vertex_idx.end());
    }
  }

  if (loop_container.added_any) {
    for (auto& loop : violating_loops) {
      std::vector<size_t> idxs;
      idxs.reserve(loop.size());
      std::transform(
          loop.begin(),
          loop.end(),
          std::back_inserter(idxs),
          [&](par::RRP::Edge& ori_idx) { return ori_idx.vertex_idx; });
      loop_container.loops.erase(idxs);
    }
    if(loop_container.loops.empty()){
      raise(SIGTRAP);
    }
    violating_loops.reserve(violating_loops.size()
                            + loop_container.loops.size());
    for (const auto& loop : loop_container.loops) {
      violating_loops.emplace_back();
      violating_loops.back().reserve(loop.size());
      auto find_edge=[&in](int prev_idx, int next_idx){
        for (auto& edge : in.vertices[prev_idx].fanouts) {
          if (edge.vertex_idx == next_idx) {
            return edge;
          }
        }
        raise(SIGTRAP);
        return par::RRP::Edge();
      };
      //First vertex
      violating_loops.back().emplace_back(find_edge(loop.back(), loop.front()));
      //Rest of the vertices
      for (size_t i = 1; i < loop.size(); i++) {
        violating_loops.back().emplace_back(find_edge(loop[i - 1], loop[i]));
      }
    }
  }

  return loop_container.added_any;
}