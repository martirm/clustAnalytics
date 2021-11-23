#include <Rcpp.h>
#include <map>
#include <vector>
#include "graph.h"
using namespace Rcpp;


//[[Rcpp::export]]
NumericVector TPR_coms_Rcpp(const IntegerVector& triangles,
                            const IntegerVector& com){
    // triangles contains the indices of the triangle vertices in succession
    // (the first three vertices belong to the first triangle, etc.)
    // com contains the communtiy index of each vertex
    // n is the number of vertices

    int n_coms = Rcpp::max(com);
    int n = com.size(), triangles_length = triangles.size();
    std::vector<bool> has_triangles (n, false);
    for (int i=0; i<triangles_length; ++i){
        int a = triangles[i++]-1, b = triangles[i++]-1, c = triangles[i]-1;
        if (com[a] == com[b] and com[a] == com[c]){
            has_triangles[a] = true;
            has_triangles[b] = true;
            has_triangles[c] = true;
        }
    }
    NumericVector TPR_coms(n_coms, 0.0), com_size(n_coms, 0.0);
    for (int i=0; i<n; ++i){
        int c = com[i] - 1; // from 1 to 0-based indices
        if (has_triangles[i]) TPR_coms[c]++;
        com_size[c]++;
    }
    //Rcout << TPR_coms << std::endl;
    //Rcout << com_size << std::endl;
    return TPR_coms / com_size;
}







//not sure this works correctly. For now use igraph implementation instead.
double triangle_participation_ratio(const Graph& g){
    int n = g.get_order();
    std::vector<bool> has_triangle(n, false); //for each vertex, indicate whether it participates in a triangle
    for (int v=0; v<n; ++v){
        if (has_triangle[v]) continue;

        auto adj = g.vertex_adjacencies(v);
        for (auto e1=adj.begin(); e1!=adj.end(); ++e1){
            int a = e1->first;
            for (auto e2=std::next(e1); e2!=adj.end(); ++e2){
                int b = e2->first;
                if (g.adjacent(a, b)){
                    has_triangle[v] = true;
                    has_triangle[a] = true;
                    has_triangle[b] = true;
                    continue;
                }
            }
        }
    }
    double ratio = 0;
    for (int i=0; i<n; ++i){
        //Rcout << (1 ? has_triangle[i] : 0) << " ";
        if (has_triangle[i]) ++ratio;
    }
    //Rcout << std::endl;
    ratio /= n;
    return ratio;
}

//[[Rcpp::export]]
double triangle_participation_ratio_Rcpp(const NumericMatrix& EdgeList){
    Graph g = Graph_from_edge_list(EdgeList);
    return triangle_participation_ratio(g);
}
