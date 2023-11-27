#include <cstdio>
#include <cstddef>
#include <vector>
namespace par{
//Correspond to single output of a cluster (a multi-output cluster /cell correspond to multiple vertices)
namespace RRP {
    struct Hgraph_Id{
        size_t cluster_id;
    };
    struct Fanin{
        unsigned int reg_cnt;
        unsigned int comb_delay;
        size_t driver_id;
    };
    struct Edge{
        size_t vertex_idx;
        size_t fanin_idx;
    };
    struct Vertex{
        //correspond to id in hgraph ds returned by coarsensing in tritonpart, opaque field
        Hgraph_Id foreign_id;
        size_t id;
        std::vector<Fanin> fanins;
        std::vector<Edge> fanouts;
        size_t edge_id_start_idx;

    };
    struct Graph{
        std::vector<Vertex> vertices;
        std::vector<size_t> primary_input_id;
        std::vector<size_t> primary_output_id;
        size_t max_component_cnt;
        int fixed_crossing_delay;
        int proportional_crossing_delay;
        int ifc_size;
        int max_delay;
        int max_crossing;
        int num_partitions;
        int total_num_edges;
    };
    void populate_fanout(RRP::Graph& g);
    Graph load_delay_graph_from_txt(const char* file_name);
    std::vector<std::vector<size_t>> load_clusters_from_txt(const char* file_name);
    void print_graph(RRP::Graph& g,FILE* ostream);
};

};