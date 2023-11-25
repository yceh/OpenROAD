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

    };
    struct Graph{
        std::vector<Vertex> vertices;
        std::vector<size_t> primary_input_id;
        std::vector<size_t> primary_output_id;
        size_t max_component_cnt;
    };
};

};