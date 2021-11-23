#include <Rcpp.h>
#include <map>
#include "graph.h"
using namespace Rcpp;


double local_clustering_coefficient(const Graph& g, const int vertex){
    // this assumes positive weights for the edges, and that weight 0 is
    // equivalent to not having an edge
    std::multimap<double,bool> sorted_triangles;
    // the key indicates at which weight a triangle or a connected vertex triple
    // appears. This is the maximum of its edge weights.
    // The value will be true for triangles, and false for triplets.
    std::map<int,double> adjacencies = g.vertex_adjacencies(vertex);
    for (auto i=adjacencies.begin(); i!=adjacencies.end(); ++i){
        int a = i->first;
        double w_va = i->second;
        for (auto j = std::next(i); j!=adjacencies.end(); ++j){
            int b = j->first;
            double w_vb = j->second;
            double max_vavb = std::max(w_va, w_vb);
            sorted_triangles.insert({max_vavb, false});

            //now find out if the vertex triple forms a closed triangle
            double w_ab = g.get_weight(a, b);
            if (w_ab > 0){
                sorted_triangles.insert({std::max(max_vavb, w_ab), true});
            }
        }
    }

    // Now we just need to integrate
    int n_triangles = 0;
    int n_triples = 0;
    double sum = 0;
    double current_weight = g.get_upper_bound(), new_weight;
    for (auto it=sorted_triangles.rbegin(); it!=sorted_triangles.rend(); ++it){
        new_weight = it->first;
        if (current_weight != new_weight and n_triples>0){
            sum += (current_weight-new_weight) * (n_triangles/double(n_triples));
            current_weight = new_weight;
        }
        if (it->second == true) n_triangles++;
        else n_triples++;
    }
    if (n_triples > 0){

        sum += current_weight * (n_triangles/double(n_triples)); //last interval of the integral (from 0 to smallest weight)
    }
    return sum/g.get_upper_bound();
}


double clustering_coefficient(const Graph& g){
    double cc = 0;
    int n = g.get_order();
    for (int i=0; i<n; ++i){
        cc += local_clustering_coefficient(g, i);
    }
    return cc/n;
}

double transitivity(const Graph& g){
    // this assumes positive weights for the edges, and that weight 0 is
    // equivalent to not having an edge
    std::multimap<double,bool> sorted_triangles;
    // the key indicates at which weight a triangle or a connected vertex triple
    // appears. This is the maximum of its edge weights.
    // The value will be true for triangles, and false for triplets.
    int n = g.get_order();
    for (int v=0; v<n; ++v){
        // we look for connected triples centered around v
        std::map<int,double> adjacencies = g.vertex_adjacencies(v);
        for (auto i=adjacencies.begin(); i!=adjacencies.end(); ++i){
            // the upper_bound function starts the iteratror at vertices with
            // index >v to prevent duplicates
            int a = i->first;
            double w_va = i->second;
            for (auto j = std::next(i); j!=adjacencies.end(); ++j){
                int b = j->first;
                double w_vb = j->second;
                double max_vavb = std::max(w_va, w_vb);
                sorted_triangles.insert({max_vavb, false}); //insert AVB triple

                //now find out if the vertex triple forms a closed triangle
                double w_ab = g.get_weight(a, b);
                if (w_ab > 0){
                    sorted_triangles.insert({std::max(max_vavb, w_ab), true}); //insert triangle
                    //sorted_triangles.insert({std::max(w_va, w_ab), false}); //insert VAB triple
                    //sorted_triangles.insert({std::max(w_vb, w_ab), false}); //insert VBA triple
                }
            }
        }
    }

    // Now we just need to integrate
    int n_triangles = 0;
    int n_triples = 0;
    double sum = 0;
    double current_weight = g.get_upper_bound(), new_weight;
    for (auto it=sorted_triangles.rbegin(); it!=sorted_triangles.rend(); ++it){
        new_weight = it->first;
        if (current_weight != new_weight and n_triples>0){
            sum += (current_weight-new_weight) * (n_triangles / double(n_triples));
            current_weight = new_weight;
        }
        if (it->second == true) n_triangles++;
        else n_triples++;
    }
    if (n_triples > 0){
        sum += current_weight * (n_triangles / double(n_triples)); //last interval of the integral (from 0 to smallest weight)
    }

    return sum/g.get_upper_bound();
}

//[[Rcpp::export]]
double clustering_coefficient_Rcpp(const NumericMatrix& EdgeList,
                              double lower_bound=0, double upper_bound=1){
    Graph g = Graph_from_edge_list(EdgeList, lower_bound, upper_bound);
    return (clustering_coefficient(g));
}

//[[Rcpp::export]]
double transitivity_Rcpp(const NumericMatrix& EdgeList,
                              double lower_bound=0, double upper_bound=1){
    Graph g = Graph_from_edge_list(EdgeList, lower_bound, upper_bound);
    return (transitivity(g));
}
