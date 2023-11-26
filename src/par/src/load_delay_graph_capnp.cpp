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
