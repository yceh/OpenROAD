#include <vector>
#include "Utilities.h"
#include "graph.hpp"
namespace RRP {
    // clusters: index of vertexes in each cluster
    // Abstract each output of a cluster into a single vertex 
    // and compute the edges after abstraction using algorithm WD 
    void make_abstracted_graph(const Graph& in,  const std::vector<std::vector<size_t>>& clusters,Graph& out);

    //For initial partitioning:
    //Provide initial sets of loops (by partitioning vetexes into different loops)
    void tarjan_loops(const Graph& in,std::vector<std::vector<Edge>>& loops);

    //The partitions each cluster is in
    struct Each_Cluster_Solution{
        std::vector<int> partitions;
        
    };
    typedef std::vector<Each_Cluster_Solution> Partition_Solution_t;
    // Finding loops with timing violations that should be added to the formulation using FEAS, and append to loops 
    void FEAS_violating_loops(const Graph& in,const Partition_Solution_t& sol_in, std::vector<std::vector<Edge>>& loops);

    //Call OR-Tools to partition
    void ilp_part(const Graph& in,Partition_Solution_t& sol_out, const std::vector<std::vector<Edge>>& loops,const Matrix<float>& upper_block_balance,
    const par::Matrix<float>& lower_block_balance);

}