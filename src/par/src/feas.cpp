#include <algorithm>
#include <boost/container_hash/extensions.hpp>
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <unordered_set>
#include <vector>

#include "delay_graph.capnp.h"
#include "replicate_retime_partition.hpp"

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
  std::vector<par::RRP::Vertex> vertices;
  std::vector<int> reg_cnt_after_retime;
};
void construct_graph_after_replication(
    const par::RRP::Partition_Solution_t& partition_solution,
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
    void add_vertex(size_t ori_idx, int partition_id)
    {
      auto idx = ori_idx * num_partitions + partition_id;
      if (new_idxes[idx] == SIZE_MAX) {
        new_idxes[idx] = g_out.vertices.size();
        g_out.vertices.emplace_back();
        g_out.attributes.emplace_back(ori_idx, partition_id);
        g_out.vertices.back().id = new_idxes[idx];
      }
    }
  };
  Idx_Map_t idx_map(in, g_out);
  g_out.vertices.reserve(in.num_partitions * in.vertices.size());
  g_out.attributes.reserve(in.num_partitions * in.vertices.size());
  for (const auto& edge : partition_solution.each_edge_solution) {
    for (const auto& partition_crossings : edge.partitions) {
      idx_map.add_vertex(edge.src_idx, partition_crossings.first);
      idx_map.add_vertex(edge.dst_idx, partition_crossings.second);
    }
  }
  for (const auto& edge : partition_solution.each_edge_solution) {
    auto reg_cnt = UINT16_MAX;
    auto comb_delay = UINT16_MAX;
    for (const auto& fanin : in.vertices[edge.dst_idx].fanins) {
      if (fanin.driver_id == edge.src_idx) {
        reg_cnt = fanin.reg_cnt;
        comb_delay = fanin.comb_delay;
        break;
      }
    }
    if (reg_cnt == UINT16_MAX || comb_delay == UINT16_MAX) {
      raise(SIGTRAP);
    }

    for (const auto& partition_crossings : edge.partitions) {
      auto src_idx = idx_map(edge.src_idx, partition_crossings.first);
      auto dst_idx = idx_map(edge.dst_idx, partition_crossings.second);
      g_out.vertices[src_idx].fanouts.emplace_back();
      g_out.vertices[src_idx].fanouts.back().vertex_idx = dst_idx;
      g_out.vertices[src_idx].fanouts.back().fanin_idx
          = g_out.vertices[dst_idx].fanins.size();
      g_out.vertices[dst_idx].fanins.emplace_back();

      auto new_comb_delay = comb_delay;
      if (partition_crossings.first != partition_crossings.second) {
        new_comb_delay
            += in.fixed_crossing_delay
               + in.proportional_crossing_delay
                     * partition_solution
                           .crossing_count[partition_crossings.first]
                                          [partition_crossings.second];
      }
      // If the comb delay containt multiple clock cycles, subtract them from
      // the register count (it has to if timing is satisfied)
      auto reg_cnt_substracted
          = std::min(reg_cnt, new_comb_delay / partition_solution.clock_period);
      new_comb_delay -= reg_cnt_substracted * partition_solution.clock_period;
      auto new_reg_cnt = reg_cnt - reg_cnt_substracted;
      g_out.vertices[dst_idx].fanins.back().reg_cnt = new_reg_cnt;
      g_out.vertices[dst_idx].fanins.back().comb_delay = new_comb_delay;
      g_out.vertices[dst_idx].fanins.back().driver_id = src_idx;
    }
  }

  size_t edge_idx = 0;
  for (size_t idx = 0; idx < g_out.vertices.size(); ++idx) {
    g_out.vertices[idx].edge_id_start_idx = edge_idx;
    edge_idx += g_out.vertices[idx].fanins.size();
  }
  g_out.reg_cnt_after_retime.resize(edge_idx);
}

// Topological traversal
struct Out_Degree_Compute
{
  Retiming_Graph& g;
  Out_Degree_Compute(Retiming_Graph& g) : g(g) {}
  size_t operator()(size_t idx) const
  {
    size_t fanout_cnt = 0;
    auto src_time_offset = g.attributes[idx].time_offset;
    for (const auto& fanout : g.vertices[idx].fanouts) {
      auto dst_vertex_idx = fanout.vertex_idx;
      const auto& dst_vertex = g.vertices[dst_vertex_idx];
      const auto& dst_attr = g.attributes[dst_vertex_idx];
      auto dst_time_offset = dst_attr.time_offset;
      auto& edge = dst_vertex.fanins[fanout.fanin_idx];
      auto& new_reg_cnt = g.reg_cnt_after_retime[dst_vertex.edge_id_start_idx
                                                 + fanout.fanin_idx];
      new_reg_cnt = (float) edge.reg_cnt + dst_time_offset - src_time_offset;
      if (new_reg_cnt < 0) {
        raise(SIGTRAP);
      }
      if (new_reg_cnt > 0) {
        ++fanout_cnt;
      }
    }
    return fanout_cnt;
  }
};

struct Time_Shifter
{
  Retiming_Graph& g;
  void operator()(size_t idx) { g.attributes[idx].time_offset++; }
};
typedef std::vector<par::RRP::Fanin> Path_t;
// Make smallest driver id first
void cannonicalize_loop(Path_t& path)
{
  auto min_iter = std::min_element(
      path.begin(), path.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.driver_id < rhs.driver_id;
      });
  std::rotate(path.begin(), min_iter, path.end());
}
struct Time_Shift_Tracker
{
  Retiming_Graph& g;
  struct Each_Vertex_Paths
  {
    std::vector<Path_t> paths;
  };
  std::vector<Each_Vertex_Paths> each_vertex_paths;
  std::unordered_set<Path_t, boost::hash<Path_t>> loops;
  Time_Shift_Tracker(Retiming_Graph& g, bool enable)
      : g(g), each_vertex_paths(enable ? g.vertices.size() : 0)
  {
  }
  void operator()(const par::RRP::Fanin& from, size_t to)
  {
    if (!each_vertex_paths.empty()) {
      return;
    }
    for (auto path : each_vertex_paths[from.driver_id].paths) {
      auto iter = std::find(path.begin(), path.end(), from);
      if (iter != path.end()) {
        Path_t loop(iter, path.end());
        cannonicalize_loop(loop);
        loops.emplace(loop);
        path.erase(iter + 1, path.end());
      } else {
        path.push_back(from);
      }
      each_vertex_paths[to].paths.push_back(path);
    }
  }
};
struct Dec_Fanout_Op
{
  std::vector<unsigned int>& fanout_cnt;
  std::vector<size_t>& next_to_traverse;
  void operator()(size_t idx)
  {
    if (--fanout_cnt[idx]) {
      next_to_traverse.push_back(idx);
    }
  }
};
struct Compute_Comb_Delay
{
  Retiming_Graph& g;
  Time_Shifter& time_shifter;
  Dec_Fanout_Op& dec_fanout_op;
  std::vector<unsigned int>& comb_delay;
  int clock_period;
  Time_Shift_Tracker& tracker;
  void operator()(size_t idx)
  {
    auto& attr = g.attributes[idx];
    // auto& vertex=g.vertices[idx];
    // Iterate over all fanins, find the critical path delay
    // if the fanin edge has register, check the delay of the driver does not
    // exceed the clock period proportion
    const auto& dst_vertex = g.vertices[idx];
    auto& this_comb_delay = comb_delay[idx];
    for (size_t fanin_idx = 0; fanin_idx < dst_vertex.fanins.size();
         fanin_idx++) {
      const auto& fanin = dst_vertex.fanins[fanin_idx];
      auto& new_reg_cnt
          = g.reg_cnt_after_retime[dst_vertex.edge_id_start_idx + fanin_idx];
      unsigned int this_fanin_delay = 0;
      // No register
      if (new_reg_cnt == 0) {
        this_fanin_delay = fanin.comb_delay + comb_delay[fanin.driver_id];
        dec_fanout_op(fanin.driver_id);
      } else {
        // Check if the immediate edge have delay greater than clock period
        // borrow the remaining delay from this segement
        int shift = fanin.comb_delay - clock_period * new_reg_cnt;
        if (shift > 0) {
          this_fanin_delay = shift;
        }
      }
      if (this_fanin_delay > clock_period) {
        tracker(fanin, idx);
      }
      this_comb_delay = std::max(this_comb_delay, this_fanin_delay);
    }
    // Delay this node if fails timing
    if (this_comb_delay > clock_period) {
      time_shifter(idx);
    }
  }
};

void topologicalSort(Retiming_Graph& g,
                     int clock_period,
                     Time_Shift_Tracker& tracker)
{
  // Compute out degree
  std::vector<unsigned int> fanout_cnt(g.vertices.size());
  std::vector<size_t> next_to_traverse;
  for (size_t idx = 0; idx < g.vertices.size(); ++idx) {
    fanout_cnt[idx] = g.vertices[idx].fanouts.size();
    if (fanout_cnt[idx] == 0) {
      next_to_traverse.push_back(idx);
    }
  }
  // Iterate over all vertices
  std::vector<unsigned int> comb_delay(g.vertices.size());
  Time_Shifter time_shifter{g};
  Dec_Fanout_Op dec_fanout_op{fanout_cnt, next_to_traverse};
  Compute_Comb_Delay compute_comb_delay{
      g, time_shifter, dec_fanout_op, comb_delay, clock_period, tracker};
  while (!next_to_traverse.empty()) {
    auto idx = next_to_traverse.back();
    next_to_traverse.pop_back();
    compute_comb_delay(idx);
  }
}
par::RRP::Edge find_edge(const par::RRP::Graph& in,
                         size_t from, size_t to)
{
    for(size_t fanin_idx=0;fanin_idx<in.vertices[to].fanins.size();fanin_idx++){
        const auto& fanin=in.vertices[to].fanins[fanin_idx];
        if(fanin.driver_id==from){
            return par::RRP::Edge{to,fanin_idx};
        }
    }
    raise(SIGTRAP);
    return par::RRP::Edge{};
}
bool par::RRP::FEAS_violating_loops(
    const par::RRP::Graph& in,
    const par::RRP::Partition_Solution_t& sol_in,
    std::vector<std::vector<par::RRP::Edge>>& loops)
{
  Retiming_Graph g;
  construct_graph_after_replication(sol_in, in, g);
  {
    Time_Shift_Tracker tracker{g, false};
    for (int iter = 0; iter < g.vertices.size(); iter++) {
      topologicalSort(g, sol_in.clock_period, tracker);
    }
  }
  {
    Time_Shift_Tracker tracker{g, true};
    for (int iter = 0; iter < g.vertices.size(); iter++) {
      topologicalSort(g, sol_in.clock_period, tracker);
    }
    typedef std::vector<size_t> Loop_t;
    struct Loop_Container
    {
      std::unordered_set<Loop_t, boost::hash<Loop_t>> loops;
      bool added_any=false;
      void operator()(Loop_t::iterator begin, Loop_t::iterator end)
      {
        Loop_t loop(begin, end);
        if (loop.size() < 2) {
            raise(SIGTRAP);
        }
        auto min_element_iter = std::min_element(begin, end);
        std::rotate(loop.begin(), min_element_iter, loop.end());
        added_any|=loops.emplace(std::move(loop)).second;
      }
    };
    Loop_Container loop_container;
    loop_container.loops.reserve(loops.size() * 2);
    for(auto & ori_loops:loops){
        Loop_t loop (ori_loops.size());
        for(size_t idx=0;idx<ori_loops.size();idx++){
            loop[idx]=ori_loops[idx].vertex_idx;
        }
        loop_container.loops.insert(loop);
    }
    loop_container.added_any=false;
    for (const auto& loop : tracker.loops) {
      int total_comb_delay = 0;
      int total_reg_cnt = 0;
      for (const auto& fanin : loop) {
        total_comb_delay += fanin.comb_delay;
        total_reg_cnt += fanin.reg_cnt;
      }
      // Filter false positive due to integer reg requirement
      if (total_comb_delay > sol_in.clock_period * total_reg_cnt) {
        // Decompose loop
        // Translate to original vertex id
        Loop_t ori_vertex_idx(loop.size());
        for (size_t i = 0; i < loop.size(); i++) {
          ori_vertex_idx[i] = g.attributes[loop[i].driver_id].ori_vertex_idx;
        }
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
    }
    if (loop_container.added_any) {
        loops.clear();
        loops.reserve(loop_container.loops.size());
        for(auto & loop:loop_container.loops){
            std::vector<par::RRP::Edge> out_loop(loop.size());
            out_loop[0]=find_edge(in,loop.back(),loop.front()); 
            for(size_t idx=1;idx<loop.size();idx++){
                out_loop[idx]=find_edge(in,loop[idx-1],loop[idx]);
            }
            loops.emplace_back(std::move(out_loop));
        }
        return true;
    }
  }
    return false;
}
