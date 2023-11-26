#include "../graph.hpp"
int main(int argc, char** argv) {
    auto g=par::RRP::load_delay_graph_from_txt(argv[1]);
    par::RRP::print_graph(g,stdout);
  return 0;
}