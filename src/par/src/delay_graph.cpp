
#include <cstdio>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include "graph.hpp"
#include <fstream>
#include <vector>
using namespace par;
void RRP::populate_fanout(RRP::Graph& g)
{
  std::size_t edge_start_idx = 0;
  for (size_t vertex_idx = 0; vertex_idx < g.vertices.size(); vertex_idx++) {
    auto& vertex = g.vertices[vertex_idx];
    vertex.edge_id_start_idx=edge_start_idx;
    for (size_t fanin_idx = 0; fanin_idx < vertex.fanins.size(); fanin_idx++) {
      const auto& fanin_edge = vertex.fanins[fanin_idx];
      g.vertices[fanin_edge.driver_id].fanouts.push_back(
          RRP::Edge{vertex_idx, fanin_idx});
    }
    edge_start_idx += vertex.fanins.size();
  }
}
RRP::Graph RRP::load_delay_graph_from_txt(const char* file_name){
    std::string line;
    std::ifstream myfile (file_name);
    RRP::Graph g;
    if (myfile.is_open())
    {
        int line_cnt=0;
        while ( getline (myfile,line) )
        {
          std::stringstream ss(line);
          line_cnt++;
          std::vector<int> fields{std::istream_iterator<int>(ss),std::istream_iterator<int>()};
          if(fields.size()==2){
            //Component ID, number of vertices with that component ID
            // vertex ID assigned consecutively
            for(int i=0;i<fields[1];i++){
              RRP::Vertex v;
              v.foreign_id.cluster_id=fields[0];
              v.id=g.vertices.size();
              g.vertices.push_back(v);
            }
          }else if (fields.size()==4){
            //destination ID, driver ID, reg_cnt, comb_delay
            RRP::Fanin fanin;
            if(fields[0]>=g.vertices.size()){
              fprintf(stderr, "@line %d: dst vertex %d does not exist, only %zu vertices\n",line_cnt,fields[0],g.vertices.size());
              continue;
            }
            if(fields[1]>=g.vertices.size()){
              fprintf(stderr, "@line %d: driver vertex %d does not exist, only %zu vertices\n",line_cnt,fields[1],g.vertices.size());
              continue;
            }
            fanin.driver_id=fields[1];
            fanin.reg_cnt=fields[2];
            fanin.comb_delay=fields[3];
            g.vertices[fields[0]].fanins.push_back(fanin);
          }else{
            std::cout<<"@ line:"<<line_cnt<<"Invalid line in delay graph file:"<<line<<std::endl;
          }
          
        }
        myfile.close();
    }
    else std::cout << "Unable to open file";
    populate_fanout(g);
    return g;
}

std::vector<std::vector<size_t>> RRP::load_clusters_from_txt(const char* file_name)
{
    std::string line;
    std::ifstream myfile (file_name);
    std::vector<std::vector<size_t>> clusters;
    if (myfile.is_open())
    {
        int line_cnt=0;
        while ( getline (myfile,line) )
        {
          std::stringstream ss(line);
          line_cnt++;
          std::vector<int> fields{std::istream_iterator<int>(ss),std::istream_iterator<int>()};
          // cluster 0: vertex0, vertex1....
          std::vector<size_t> cluster;
          for (int v : fields) {
            cluster.push_back(v);
          }
          clusters.push_back(cluster);          
        }
        myfile.close();
    }
    else std::cout << "Unable to open file";

    return clusters;
}

void RRP::print_graph(RRP::Graph& g,FILE* ostream){
  fprintf(ostream,"Graph has %lu vertices\n",g.vertices.size());
  for(auto& v:g.vertices){
    fprintf(ostream,"Vertex %lu: Component id %lu, edge_start_idx %lu \n",v.id,v.foreign_id.cluster_id,v.edge_id_start_idx);
    fprintf(ostream,"Fanin: \n");
    int idx=0;
    for(auto& fanin:v.fanins){
      fprintf(ostream,"\t %d:driver %lu, reg_cnt %u, comb_delay %u\n",idx++,fanin.driver_id,fanin.reg_cnt,fanin.comb_delay);
    }
    idx=0;
    fprintf(ostream,"Fanout: \n");
    for(auto& fanout:v.fanouts){
      fprintf(ostream,"\t %d:dst %lu, fanin_idx %lu\n",idx++,fanout.vertex_idx,fanout.fanin_idx);
    }
  }
}