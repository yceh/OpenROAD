#include "replicate_retime_partition.hpp"
#include <cstddef>
#include <stack>
#include <queue>
#include <unordered_map>
#include <vector>
namespace par{
namespace RRP {

    #define undefined_index std::numeric_limits<size_t>::max()

    // create a mapping table to record the cluster number of each vertex
    void cluster_num_of_vertex(const std::vector<std::vector<size_t>>& clusters,  std::vector<size_t>& vertices_cluster_num)
    {
        size_t cluster_idx = 0;

        for (const std::vector<size_t>& cluster : clusters) {
            for (const size_t& v : cluster) {
                vertices_cluster_num[v] = cluster_idx;
            }
            cluster_idx++;
        }
    }


    struct src_inf {
        size_t vertex_idx;
        std::vector<size_t> fanins;
    };


    void Find_src_dst(const Graph& in, const size_t cluster_num, const std::vector<size_t>& cluster, const std::vector<size_t>& vertices_cluster_num, 
                        std::vector<src_inf> &src, std::vector<size_t> &dst)
    {
        for (const size_t& v : cluster) {

            // searching dst
            for (const Edge& e : in.vertices[v].fanouts) {
                const size_t& next_v = e.vertex_idx;
                if (vertices_cluster_num[next_v] != cluster_num) {
                    dst.push_back(v);
                }
            }

            // searching src
            size_t fanin_idx = 0;
            std::vector<size_t> src_fanins;
            for (const Fanin& fanin : in.vertices[v].fanins) {
                const size_t& pre_v = fanin.driver_id;
                if (vertices_cluster_num[pre_v] != cluster_num) {
                    src_fanins.push_back(fanin_idx);
                }
                fanin_idx++;
            }
            if (!src_fanins.empty()) {
                src.push_back({v, src_fanins});
            }
        }
    }

    const size_t INF = std::numeric_limits<unsigned int>::max(); // Infinity value for initialization

    struct reg_delay {
        unsigned int reg_cnt;
        int neg_delay;

    };


    // searching minimum distance from soure to each destination (smaller reg and smaller neg_delay)
    void Dijkstra (const Graph& in, const size_t cluster_num, const std::vector<size_t>& cluster, const std::vector<size_t>& vertices_cluster_num, 
                    const size_t &source, const std::vector<size_t> &dsts, std::vector<Fanin> &new_fanin)
    {
        std::unordered_map<size_t, reg_delay> distances;
        distances.reserve(cluster.size()); // Initialize distances to infinity

        // Custom comparison for the priority queue
        auto cmp = [](std::pair<size_t, reg_delay> left, std::pair<size_t, reg_delay> right) {
            if (left.second.reg_cnt != right.second.reg_cnt) {
                return left.second.reg_cnt < right.second.reg_cnt;
            } else {
                return left.second.neg_delay < right.second.neg_delay;
            }
        };

        // first : vertex_idx; second : reg_delay
        std::priority_queue<std::pair<size_t, reg_delay>, std::vector<std::pair<size_t, reg_delay>>, decltype(cmp)> pq(cmp);

        distances[source].reg_cnt = 0; // Distance from source to itself is 0
        distances[source].neg_delay = 0;
        pq.push({source, distances[source]}); // Push the source vertex and its distance

        while (!pq.empty()) {
            size_t u = pq.top().first;  
            auto u_weight=pq.top().second;        
            pq.pop();

            // Iterate through the adjacency list of the current vertex
            for (const Edge& edge : in.vertices[u].fanouts) {
                size_t next_v = edge.vertex_idx;
                // only consider the vertex inside the cluster
                if (vertices_cluster_num[next_v] == cluster_num) {
                    size_t fanin_idx = edge.fanin_idx;
                    unsigned int reg_cnt = in.vertices[next_v].fanins[fanin_idx].reg_cnt;
                    int neg_delay = -1*(in.vertices[next_v].fanins[fanin_idx].comb_delay);
                    reg_delay v_weight{u_weight.reg_cnt+reg_cnt, u_weight.neg_delay+neg_delay};
                    auto result = distances.emplace(next_v, v_weight);
                    if(result.second){
                        pq.push({next_v, v_weight});
                    }else{
                        auto & inserted_delay=result.first->second;
                        if (inserted_delay.reg_cnt>v_weight.reg_cnt || 
                                (inserted_delay.reg_cnt == v_weight.reg_cnt && inserted_delay.neg_delay > v_weight.neg_delay) ) {
                            inserted_delay = v_weight;
                            pq.push({next_v, v_weight});
                        }
                    }
                }
            }
        }

        // only return the distance regarding dsts
        for (const size_t dst : dsts) {
            std::unordered_map<size_t, reg_delay>::const_iterator got = distances.find (dst);
            Fanin n_fanin;
            if (got == distances.end() ) {
                // src does not connect to this dst
                n_fanin.comb_delay = 0;
                n_fanin.reg_cnt = INF;
                n_fanin.driver_id = source;
            } else {
                // src connects to this dst
                n_fanin.comb_delay = -1*got->second.neg_delay;
                n_fanin.reg_cnt = got->second.reg_cnt;
                n_fanin.driver_id = source;
            }
            new_fanin.push_back(n_fanin);
        }

    }

    void add_fanin_inf(const std::vector<Fanin> &new_fanin, const Fanin &src_fanin, std::vector<Fanin> &result_fanin) 
    {
        for (const Fanin n_fanin : new_fanin) {
            if (n_fanin.comb_delay == 0 && n_fanin.reg_cnt == INF) {
                result_fanin.push_back( {n_fanin.reg_cnt, 
                                            n_fanin.comb_delay, 
                                            src_fanin.driver_id}); 
            } else {
                result_fanin.push_back( {n_fanin.reg_cnt+src_fanin.reg_cnt, 
                                            n_fanin.comb_delay+src_fanin.comb_delay, 
                                            src_fanin.driver_id});                
            }
        }
    }


    void WD_refine(const Graph& in, const size_t cluster_num, const std::vector<size_t>& cluster, const std::vector<size_t> &vertices_cluster_num, 
                    const std::vector<src_inf> &srcs, const std::vector<size_t> &dsts, std::vector<std::vector<Fanin>> &dsts_new_fanins) 
    {
        std::vector<std::vector<Fanin>> srcs_based_new_fanins;
                
        for (const src_inf& src : srcs) {
            size_t src_vertex = src.vertex_idx;
            std::vector<Fanin> new_fanin;
            
            Dijkstra(in, cluster_num, cluster, vertices_cluster_num, src_vertex, dsts, new_fanin);

            // for debug----------------------------------
            // printf("src id:%lu\n", src.vertex_idx);
            // for (Fanin f : new_fanin) {
            //     printf("fanin reg cnt: %d, comb delay: %d, driver id: %lu\n", f.reg_cnt, f.comb_delay, f.driver_id);
            // }

            //--------------------------------------------


            for (const size_t fanin_idx : src.fanins) {
                Fanin src_fanin = in.vertices[src_vertex].fanins[fanin_idx];
                std::vector<Fanin> result_fanin;
                add_fanin_inf(new_fanin, src_fanin, result_fanin);
                srcs_based_new_fanins.push_back(result_fanin);
            }
        }
        for (size_t dst_idx = 0; dst_idx < dsts.size(); dst_idx++) {
            std::vector<Fanin> dst_fanins;
            for (size_t src_idx = 0; src_idx < srcs.size(); src_idx++) {
                Fanin dst_fanin = srcs_based_new_fanins[src_idx][dst_idx];
                // valid connection
                if (!(dst_fanin.reg_cnt == INF &&dst_fanin.comb_delay == 0)) {
                    dst_fanins.push_back(dst_fanin);
                }
            }
            dsts_new_fanins.push_back(dst_fanins);
        }

    }

    // remove vertice in the cluster and add new vertices 
    void create_new_graph(const std::vector<std::vector<size_t>> &all_dsts, const std::vector<std::vector<std::vector<Fanin>>> &all_new_fanins, 
                        std::vector<size_t>& mapping_table, Graph& new_graph)
    {
        // creating a new graph
        // 1. add all new vertex (from dsts of old graph)
        size_t new_vertex_id = 0;
        size_t cluster_id = 0;
        for (const std::vector<size_t>& dsts : all_dsts) {
            for (const size_t& dst : dsts) {
                // adding new vertex
                Vertex new_graph_vertex;
                new_graph_vertex.id = new_vertex_id;
                new_graph_vertex.foreign_id.cluster_id = cluster_id;
                new_graph.vertices.push_back(new_graph_vertex);
                // updating mapping table
                mapping_table.push_back(dst);   // index : new vertex id; value : old vertex id                
                
                new_vertex_id++;
            }
            cluster_id++;
        }
        // 2. update fanin edges in new graph
        std::unordered_map<size_t, size_t> reversed_mapping_table;
        for (size_t new_idx = 0; new_idx < mapping_table.size(); new_idx++) {
            reversed_mapping_table.emplace(mapping_table[new_idx], new_idx);
        }
        new_vertex_id = 0;
        for(const std::vector<std::vector<Fanin>>& fanins : all_new_fanins) {
            for(const std::vector<Fanin>& fanin_of_a_dst : fanins) {
                for (const Fanin& a_fanin : fanin_of_a_dst) {
                    Fanin n_fanin;          // fanin for new graph
                    n_fanin.comb_delay = a_fanin.comb_delay;
                    n_fanin.reg_cnt = a_fanin.reg_cnt;
                    auto reversed_got = reversed_mapping_table.find(a_fanin.driver_id);
                    if (reversed_got == reversed_mapping_table.end()) {
                        n_fanin.driver_id = undefined_index;
                    } else {
                        n_fanin.driver_id = reversed_got->second;                    
                    }
                    new_graph.vertices[new_vertex_id].fanins.push_back(n_fanin); 
                }

                new_vertex_id++;
            }
        }
    }


    // clusters: index of vertexes in each cluster
    // Abstract each output of a cluster into a single vertex 
    // and compute the edges after abstraction using algorithm WD 
    void make_abstracted_graph(const Graph& in,  const std::vector<std::vector<size_t>>& clusters, Graph& out, std::vector<size_t>& mapping_table)
    {
        //1. use WD algorithm obtain the  reg/delay
        //2. create a new vertex to replace
        //3. change the fanin vertex of output (dst)

        size_t numVertices = in.vertices.size();
        std::vector<size_t> vertices_cluster_num(numVertices);
        // create cluster number mapping table
        cluster_num_of_vertex(clusters, vertices_cluster_num);

        size_t cluster_idx = 0;
        std::vector<std::vector<size_t>> all_dsts;
        std::vector<std::vector<std::vector<Fanin>>> all_new_fanins;

        for (const std::vector<size_t>& cluster : clusters) {
            std::vector<src_inf> srcs;
            std::vector<size_t> dsts;
            std::vector<std::vector<Fanin>> new_fanins;
            Find_src_dst(in,  cluster_idx, cluster, vertices_cluster_num, srcs, dsts);

            // debug-------------
            // for (src_inf src: srcs) {
            //     printf("src vertex: %lu\n", src.vertex_idx);
            //     for (size_t fanin_idx : src.fanins) {
            //         printf("    fanin idx %lu\n", fanin_idx);
            //     }
            // }

            // for (size_t dst: dsts) {
            //     printf("dst vertex: %lu\n", dst);
            // }
            //-----------------


            WD_refine(in, cluster_idx, cluster, vertices_cluster_num, srcs, dsts, new_fanins);
            all_dsts.push_back(dsts);
            all_new_fanins.push_back(new_fanins);
            cluster_idx++;
        }

        create_new_graph(all_dsts, all_new_fanins, mapping_table, out);

    }

    // ------------tarjan--------------------------------------------

    // Helper struct to represent vertices in Tarjan's algorithm
    struct TarjanVertex {
        size_t index; // The index of the vertex
        size_t lowlink; // The lowlink value in Tarjan's algorithm
        bool onStack; // Flag to determine if the vertex is on the stack
    };

    // Recursive function for Tarjan's algorithm
    void tarjanDFS(size_t v, size_t& index, std::stack<Edge>& stack,
                   std::vector<TarjanVertex>& vertices, const Graph& graph,
                   std::vector<std::vector<Edge>>& loops) {
        vertices[v].index = index;
        vertices[v].lowlink = index;
        index++;
        // stack.push(v); ?
        vertices[v].onStack = true;

        // Traverse the fanouts of the current vertex
        for (const Edge& e : graph.vertices[v].fanouts) {
            size_t w = e.vertex_idx;
            stack.push(e);
            if (vertices[w].index == undefined_index) {
                tarjanDFS(w, index, stack, vertices, graph, loops);
                vertices[v].lowlink = std::min(vertices[v].lowlink, vertices[w].lowlink);
            } else if (vertices[w].onStack) {
                vertices[v].lowlink = std::min(vertices[v].lowlink, vertices[w].index);
            }
        }

        // If v is the root of a strongly connected component
        if (vertices[v].lowlink == vertices[v].index) {
            std::vector<Edge> loop;
            Edge loop_edge;
            size_t w;
            do {
                loop_edge = stack.top();
                stack.pop();
                w = loop_edge.vertex_idx;
                vertices[w].onStack = false;
                loop.push_back(loop_edge); // Add vertex to the loop
            } while (w != v);
            if (loop.size() > 1) {
                loops.push_back(loop); // Add the loop to the loops vector
            }
        }
    }

    //For initial partitioning:
    //Provide initial sets of loops (by partitioning vetexes into different loops)
    void tarjan_loops(const Graph& in, std::vector<std::vector<Edge>>& loops) {
        size_t numVertices = in.vertices.size();
        std::vector<TarjanVertex> vertices(numVertices, {undefined_index, 0, false});
        // std::stack<Edge> stack;
        size_t index = 0;

        for (size_t v = 0; v < numVertices; ++v) {
            if (vertices[v].index == undefined_index) {
                // Do we want every vertex start with empty stack?
                // Does it generate the replicated loop?
                std::stack<Edge> stack; 
                tarjanDFS(v, index, stack, vertices, in, loops);
            }
        }
    }



}
}