#include "graph.hpp"
#include "Hypergraph.h"
#include "TritonPart.h"
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
