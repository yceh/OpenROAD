#include "Hypergraph.h"
#include "TritonPart.h"
#include "graph.hpp"
#include "delay_graph.capnp.h"
#include <capnp/message.h>
#include <capnp/serialize-packed.h>
#include <fcntl.h>
#include <unistd.h>
#include <algorithm>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>
using namespace par;
void populate_fanout(RRP::Graph& g)
{
  for (size_t vertex_idx = 0; vertex_idx < g.vertices.size(); vertex_idx++) {
    const auto& vertex = g.vertices[vertex_idx];
    for (size_t fanin_idx = 0; fanin_idx < vertex.fanins.size(); fanin_idx++) {
      const auto& fanin_edge = vertex.fanins[fanin_idx];
      g.vertices[fanin_edge.driver_id].fanouts.push_back(
          RRP::Edge{vertex_idx, fanin_idx});
    }
  }
}
RRP::Graph* load_graph_from_file(const char* filename)
{
    RRP::Graph* g_ptr = new RRP::Graph;
  RRP::Graph& g=*g_ptr;
  g.max_component_cnt = 0;
  {
    ::capnp::ReaderOptions options;
    options.traversalLimitInWords = 1LL << 60;
    auto fd = open(filename, O_RDONLY);
    ::capnp::PackedFdMessageReader message(fd, options);

    Graph::Reader graph = message.getRoot<Graph>();

    g.vertices.resize(graph.getVertices().size());
    for (size_t vertex_idx = 0; vertex_idx < graph.getVertices().size();
         vertex_idx++) {
      const auto& vertex = graph.getVertices()[vertex_idx];
      auto& vertex_out = g.vertices[vertex_idx];
      vertex_out.id = vertex.getId();
      vertex_out.foreign_id.cluster_id = vertex.getComponentId();
      g.max_component_cnt
          = std::max(g.max_component_cnt, vertex_out.foreign_id.cluster_id);
      vertex_out.fanins.resize(vertex.getFanins().size());
      for (size_t fanin_idx = 0; fanin_idx < vertex.getFanins().size();
           fanin_idx++) {
        const auto& fanin = vertex.getFanins()[fanin_idx];
        auto& fanin_out = vertex_out.fanins[fanin_idx];
        fanin_out.driver_id = fanin.getDriverId();
        fanin_out.comb_delay = fanin.getCombDelay();
        fanin_out.reg_cnt = fanin.getRegCnt();
      }
    }
    close(fd);
  }

  g.max_component_cnt++;
  // Build Fanout:
  populate_fanout(g);
  return g_ptr;
}
void par::TritonPart::hyperGraph_from_delay_graph(const RRP::Graph*g_ptr){
    const RRP::Graph& g=*g_ptr;
    num_hyperedges_=g.vertices.size();
    num_vertices_=g.max_component_cnt;
    hyperedges_.clear();
    hyperedges_.resize(g.vertices.size());
    for (const auto& vertex : g.vertices) {
        auto& hyperedge_out = hyperedges_[vertex.id];
        hyperedge_out.push_back(vertex.foreign_id.cluster_id);
        for (const auto& fanout : vertex.fanouts) {
            hyperedge_out.push_back(g.vertices[fanout.vertex_idx].foreign_id.cluster_id);
        }
        std::sort(hyperedge_out.begin(), hyperedge_out.end());
        hyperedge_out.erase(std::unique(hyperedge_out.begin(), hyperedge_out.end()), hyperedge_out.end());
    }
    hyperedge_weights_=std::vector<std::vector<float>>(num_hyperedges_,std::vector<float>(hyperedge_dimensions_,1.0));
    vertex_weights_=std::vector<std::vector<float>>(num_vertices_,std::vector<float>(hyperedge_dimensions_,1.0));
    fixed_attr_.clear();
    community_attr_.clear();
    group_attr_.clear();
    placement_attr_.clear();
    original_hypergraph_ = std::make_shared<Hypergraph>(vertex_dimensions_,
                                                      hyperedge_dimensions_,
                                                      placement_dimensions_,
                                                      hyperedges_,
                                                      vertex_weights_,
                                                      hyperedge_weights_,
                                                      fixed_attr_,
                                                      community_attr_,
                                                      placement_attr_,
                                                      logger_);
}
