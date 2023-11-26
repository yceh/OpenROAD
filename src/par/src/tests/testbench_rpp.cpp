#include <cstddef>
#include "replicate_retime_partition.hpp"

using namespace par;

void build_test_graph(RRP::Graph& test_graph);

// for displaying 
void show_graph(RRP::Graph& graph, std::vector<size_t>& mapping_table);

int main()
{
    RRP::Graph test_graph, new_graph;
    std::vector<size_t> mapping_table;
    std::vector<std::vector<size_t>> clusters;
    std::vector<size_t> cluster0, cluster1;

    build_test_graph(test_graph);
    cluster0.push_back(0);
    cluster0.push_back(1);
    cluster0.push_back(3);
    cluster0.push_back(4);
    cluster1.push_back(2);
    cluster1.push_back(5);
    clusters.push_back(cluster0);
    clusters.push_back(cluster1);

    RRP::make_abstracted_graph(test_graph, clusters, new_graph, mapping_table);

    show_graph(new_graph, mapping_table);

 
    return 0;
}

void build_test_graph (RRP::Graph& test_graph) 
{
    RRP::Vertex v0, v1, v2, v3, v4, v5;
    v0.id = 0;
    v0.fanins.push_back({0, 2, 3});
    v0.fanouts.push_back({1, 0});
    test_graph.vertices.push_back(v0);
    v1.id = 1;
    v1.fanins.push_back({1, 0, 0});
    v1.fanouts.push_back({2, 0});
    v1.fanouts.push_back({4, 0});
    test_graph.vertices.push_back(v1);
    v2.id = 2;
    v2.fanins.push_back({1, 2, 1});
    v2.fanouts.push_back({5, 0});
    test_graph.vertices.push_back(v2);
    v3.id = 3;
    v3.fanins.push_back({1, 1, 4});
    v3.fanouts.push_back({0, 0});
    test_graph.vertices.push_back(v3);
    v4.id = 4;
    v4.fanins.push_back({0,1, 1});
    v4.fanins.push_back({1, 1, 5});
    v4.fanouts.push_back({3, 0});
    test_graph.vertices.push_back(v4);
    v5.id = 5;
    v5.fanins.push_back({1, 3, 2});
    v5.fanouts.push_back({4, 1});
    test_graph.vertices.push_back(v5);
}

void show_graph(RRP::Graph& graph, std::vector<size_t>& mapping_table)
{
    for (RRP::Vertex& v : graph.vertices) {
        printf("new vertex : %lu;  old vertex : %lu;\n", v.id, mapping_table[v.id]);
        for (RRP::Fanin& fanin : graph.vertices[v.id].fanins) {
            printf("fanin: reg: %d, delay: %d, driver id: %lu\n", fanin.reg_cnt, fanin.comb_delay, fanin.driver_id);
        }
        for (RRP::Edge& fanout : graph.vertices[v.id].fanouts) {
            printf("fanout: vertex_id: %lu, fanin_id: %lu\n", fanout.vertex_idx, fanout.fanin_idx);
        }
    }
}