#include <cstddef>
#include <utility>
#include <vector>
#include "Utilities.h"
#include "graph.hpp"
namespace par{
namespace RRP {
    // clusters: index of vertexes in each cluster
    // Abstract each output of a cluster into a single vertex 
    // and compute the edges after abstraction using algorithm WD 
    void make_abstracted_graph(const Graph& in,  const std::vector<std::vector<size_t>>& clusters,Graph& out, std::vector<size_t>& mapping_table);

    //For initial partitioning:
    //Provide initial sets of loops (by partitioning vetexes into different loops)
    void tarjan_loops(const Graph& in,std::vector<std::vector<Edge>>& loops);

    //The partitions each cluster is in
    struct Each_Edge_Solution{
        std::vector<std::pair<int, int>> partitions;
        size_t src_idx;
        size_t dst_idx;
    };
    struct Component_To_Output_Solution{
        std::vector<std::pair<int, int>> partitions;
    };
    struct Partition_Solution_t{
        std::vector<Component_To_Output_Solution> component_to_output_solutions;
        std::vector<Each_Edge_Solution> each_edge_solution;
        int clock_period;
        std::vector<std::vector<int>> crossing_count;
    };
    // Finding loops with timing violations that should be added to the formulation using FEAS, and append to loops 
    bool FEAS_violating_loops(const Graph& in,const Partition_Solution_t& sol_in, std::vector<std::vector<Edge>>& loops);

    //Call OR-Tools to partition
    void ilp_part(const Graph& in,Partition_Solution_t& sol_out, const std::vector<std::vector<Edge>>& loops,std::vector<size_t> max_capacity);

    void find_initial_loops(const Graph& g, std::vector<std::vector<Edge>>& loops);
}
}