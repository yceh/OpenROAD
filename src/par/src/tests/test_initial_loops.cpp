#include <cstddef>
#include "../replicate_retime_partition.hpp"

using namespace par;

void build_test_graph(RRP::Graph& test_graph);

// for displaying 
void show_graph(RRP::Graph& graph, std::vector<size_t>& mapping_table);
void show_loops(RRP::Graph& graph, std::vector<std::vector<RRP::Edge>>& loops);

int main(int argc, char** argv)
{
    RRP::Graph test_graph, new_graph;
    std::vector<size_t> mapping_table;
    std::vector<std::vector<size_t>> clusters;
    std::vector<std::vector<RRP::Edge>> loops;

    test_graph=par::RRP::load_delay_graph_from_txt(argv[1]);
    
    // test for tarjan algorithm
    RRP::find_initial_loops(test_graph, loops);
    show_loops(test_graph, loops);

 
    return 0;
}

void show_loops(RRP::Graph& graph, std::vector<std::vector<RRP::Edge>>& loops)
{
    size_t loop_id = 0;
    for (const std::vector<RRP::Edge>& loop : loops) {
        printf("loop#%lu\n", loop_id);
        for (RRP::Edge e : loop) {
            size_t driver_id = graph.vertices[e.vertex_idx].fanins[e.fanin_idx].driver_id;
            size_t dst_id = e.vertex_idx;
            printf("        vertex #%lu to vertex #%lu\n", driver_id, dst_id);
        }
        loop_id++;
    }
}