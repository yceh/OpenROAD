#include <absl/container/flat_hash_map.h>

#include <csignal>
#include <cstddef>
#include <cstdlib>
#include <vector>

#include "replicate_retime_partition.hpp"

typedef std::vector<par::RRP::Edge> Loop_t;
typedef std::vector<Loop_t> Loops_t;
namespace Find_Initial_Loop_NS {
struct Reg_Delay
{
  int reg_cnt;
  int delay;
  Reg_Delay(int reg_cnt, int delay) : reg_cnt(reg_cnt), delay(delay) {}
  Reg_Delay() = default;
  Reg_Delay(const par::RRP::Fanin& fanin)
      : reg_cnt(fanin.reg_cnt), delay(fanin.comb_delay)
  {
  }
  bool operator<(const Reg_Delay& other) const
  {
    return reg_cnt < other.reg_cnt
           || (reg_cnt == other.reg_cnt && delay > other.delay);
  }
  Reg_Delay& operator+=(const Reg_Delay& other)
  {
    reg_cnt += other.reg_cnt;
    delay += other.delay;
    return *this;
  }
  Reg_Delay operator+(const Reg_Delay& other) const
  {
    Reg_Delay tmp = *this;
    tmp += other;
    return tmp;
  }
};

// Recording which vertexes(low_link) on the call stack can be reached(form a
// loop), and the fanout to reach it with min reg_delay.
struct Low_Link_Attr
{
  Reg_Delay reg_delay;
  int fanout_min_reg_delay;
};

using Low_Links_t = absl::flat_hash_map<int, Low_Link_Attr>;
class Per_Vertex_Traversal_Attr
{
  Low_Links_t
      low_links;  // which vertex on the call stack before me can be reached
  long dfs_idx;   // 0: not visited, >0: visited,on call stack <0: visited,but
                  // not on call stack (not reachable from the current DFS root)

 public:
  static int next_dfs_idx;
  Per_Vertex_Traversal_Attr() : dfs_idx(0) {}

  // RAII for preventing to forget to mark vertex not on the stack
  struct Visiting
  {
    Per_Vertex_Traversal_Attr& self;
    bool new_visit;
    ~Visiting()
    {
      if (new_visit) {
        if (!self.is_on_call_stack()) {
          raise(SIGTRAP);
        }
        self.done_visit();
      }
    }
  };
  Visiting visit()
  {
    if (dfs_idx == 0) {
      dfs_idx = ++next_dfs_idx;
      return Visiting{*this, true};
    }
    return Visiting{*this, false};
  }

  void done_visit()
  {
    if (dfs_idx > 0) {
      dfs_idx = -dfs_idx;
    }
    if (dfs_idx >= 0) {
      raise(SIGTRAP);
    }
  }

  long get_dfs_idx() const { return abs(dfs_idx); }

  bool is_visited() const { return dfs_idx != 0; }

  bool is_on_call_stack() const { return dfs_idx > 0; }

  void low_link_reserve(int other_lowlink_cnt)
  {
    low_links.reserve(low_links.size() + other_lowlink_cnt + 1);
  }

  Low_Links_t& get_low_links() { return low_links; }

  void add_low_link(int low_link, int fanout_idx, const Reg_Delay& reg_delay)
  {
    auto res
        = low_links.emplace(low_link, Low_Link_Attr{reg_delay, fanout_idx});
    if (!res.second) {
      // There are multiple ways to reach the same vertex on the call stack.
      // Choose the one with min reg_delay
      if (res.first->second.reg_delay < reg_delay) {
        res.first->second.reg_delay = reg_delay;
        res.first->second.fanout_min_reg_delay = fanout_idx;
      }
    }
  }
};
#define DEBUG_INITIAL_LOOP
int Per_Vertex_Traversal_Attr::next_dfs_idx = 0;
struct DFS_Visitor
{
  const par::RRP::Graph& g;
  std::vector<Per_Vertex_Traversal_Attr> visit_attributes;
  Loops_t& output;
  std::vector<int> dfs_idx_to_vertex_idx_map;
  DFS_Visitor(const par::RRP::Graph& g, Loops_t& output)
      : g(g), visit_attributes(g.vertices.size()), output(output)
  ,dfs_idx_to_vertex_idx_map(g.vertices.size()+1)
  {
  }

  void accumulate_loop(size_t idx, int low_link_idx)
  {
    auto& attr = visit_attributes[idx];
    if(attr.get_dfs_idx()==low_link_idx){
      return;
    }
    auto low_links = attr.get_low_links();
    auto iter = low_links.find(low_link_idx);
    if (iter == low_links.end()) {
      raise(SIGTRAP);
    }
    auto& loop = output.back();
    size_t fanout_idx = iter->second.fanout_min_reg_delay;
    auto& fanout=g.vertices[idx].fanouts[fanout_idx];
    loop.push_back(fanout);

    accumulate_loop(fanout.vertex_idx,
                    low_link_idx);
    low_links.erase(iter);
  }
  // visit a vertex, if already visited,
  // return whether it is on the call stack before its caller,
  //  if not visited, return false
  bool operator()(int idx)
  {
    auto& attr = visit_attributes[idx];
    auto visiting = attr.visit();
    if (!visiting.new_visit) {
      return attr.is_on_call_stack();
    }
    auto this_dfs_idx=attr.get_dfs_idx();
  #ifdef DEBUG_INITIAL_LOOP
    dfs_idx_to_vertex_idx_map[this_dfs_idx]=idx;
  #endif
    auto& vertex = g.vertices[idx];
    for (size_t fanout_idx = 0; fanout_idx < vertex.fanouts.size();
         ++fanout_idx) {
      auto& fanout = vertex.fanouts[fanout_idx];
      // child parent refer to the spanning tree built by DFS
      auto child_idx = fanout.vertex_idx;
      auto reached_parent = (*this)(child_idx);
      auto& child_vertex = g.vertices[child_idx];
      auto& to_child_edge = child_vertex.fanins[fanout.fanin_idx];
      auto to_child_reg_delay = Reg_Delay(to_child_edge);
      if (reached_parent) {
        // Found back edge, the last edge of a loop
        auto& child_attr=visit_attributes[child_idx];
        attr.add_low_link(child_attr.get_dfs_idx(), fanout_idx, to_child_reg_delay);
        if (!child_attr.get_low_links().empty()) {
          raise(SIGTRAP);
        }
      }
      // Accumulate the low link of the child that are still on the call stack
      // (will form a loop in this iteration of DFS)
      auto& child_attr = visit_attributes[child_idx];
      
      for (auto& child_low_link : child_attr.get_low_links()) {
        if (visit_attributes[dfs_idx_to_vertex_idx_map[child_low_link.first]].is_on_call_stack()) {
          attr.add_low_link(
              child_low_link.first,
              fanout_idx,
              child_low_link.second.reg_delay + to_child_reg_delay);
        }
      }
      
    }

    // See whether we complete a loop (some low link = idx)
    auto& low_links = attr.get_low_links();
    auto iter = low_links.find(this_dfs_idx);
    if (iter != low_links.end()) {
      // Found a loop
      output.emplace_back();
      auto& loop = output.back();
      size_t fanout_idx = iter->second.fanout_min_reg_delay;
      auto& fanout=vertex.fanouts[fanout_idx];
      loop.push_back(fanout);
      accumulate_loop(fanout.vertex_idx, this_dfs_idx);
      low_links.erase(iter);
    }

    //DEBUG
    for(auto& low_link : low_links) {
      if (low_link.first >= this_dfs_idx) {
        raise(SIGTRAP);
      }
    }
    return false;
  }
};

};  // namespace Find_Initial_Loop_NS

void par::RRP::find_initial_loops(const par::RRP::Graph& g, std::vector<std::vector<Edge>>& loops)
{
  Find_Initial_Loop_NS::Per_Vertex_Traversal_Attr::next_dfs_idx = 0;
  Find_Initial_Loop_NS::DFS_Visitor visitor(g, loops);
  for (size_t i = 0; i < g.vertices.size(); ++i) {
    visitor(i);
  }
}