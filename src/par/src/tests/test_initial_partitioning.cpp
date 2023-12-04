#include <cstdint>
#include <cstdio>
#include "../replicate_retime_partition.hpp"
int main(int argc, char** argv)
{
  auto g = par::RRP::load_delay_graph_from_txt(argv[1]);
  for (auto& v : g.vertices) {
    v.foreign_id.cluster_id = v.id;
    for(auto& fanin:v.fanins){
      fanin.comb_delay*=10;
    }
  }
  std::vector<std::vector<par::RRP::Edge>> loops;
  par::RRP::find_initial_loops(g, loops);
  g.num_partitions = atoi(argv[2]);
  std::vector<size_t> partition_sizes(g.num_partitions,
                                      g.vertices.size() / g.num_partitions + 2);
  g.max_delay=1000000;
  g.ifc_size=2;
  g.proportional_crossing_delay=1;
  g.fixed_crossing_delay=1;
  g.max_component_cnt=g.vertices.size();
  par::RRP::Partition_Solution_t sol;
  for (int i=0;true;i++) {
    par::RRP::ilp_part(g, sol, loops, partition_sizes,std::vector<int64_t>(g.vertices.size(),1));
    fprintf(stderr, "iter %d: \n",i);
    fprintf(stderr, "\tclock period: %d\n",sol.clock_period);
    for (size_t out_vertex_idx=0; out_vertex_idx<g.vertices.size(); out_vertex_idx++) {
      fprintf(stderr, "\tvertex %zu: \n",out_vertex_idx);
      for (const auto & part: sol.component_to_output_solutions[out_vertex_idx].partitions) {
        fprintf(stderr, "\t\t%d-> %d\n",part.first,part.second);
      }
      fprintf(stderr, "\n");
    }

    if (!par::RRP::find_violating_loops(sol, g, loops)) {
      break;
    }
  }
}