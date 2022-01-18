#ifndef GRAPH_H // include guard
#define GRAPH_H

#include <map>
#include <vector>
#include <Rcpp.h>
#include "svector.h"

struct Edge{
    int a;
    int b;
    double weight;
};

//Undirected weighted graph
class Graph {
        int order;
        int size;
        double lower_bound;
        double upper_bound;
        bool directed;
        std::map<std::pair<int,int>, double> edge_list_m; //keys are pairs of ints indicating the two ends of the edge. Value is a double with its weight
        std::vector<std::map<int,double> > adjacencies_list; //for every vertex, map of its adjacencies (pairs of vertex index, and edge weight)
        SVector<std::pair<int,int> > sampling_vector; //contains all edges. Used to sample edges uniformly.

        void update_size();

    public:
        Graph(std::vector<Edge> v, int graph_order, double lb, double ub, bool d);
        Graph(std::vector<Edge> v, int graph_order, double lb, double ub) : Graph (v, graph_order, lb, ub, false) {};
        Graph(std::vector<Edge> v, int graph_order) : Graph(v, graph_order, 0, DBL_MAX) {};
        Rcpp::NumericMatrix numericmatrix_edgelist();
        int get_size() const;
        int get_order() const;
        double get_upper_bound() const;
        double get_lower_bound() const;
        double get_weight(int a, int b) const;
        bool adjacent(int a, int b) const;
        void delete_edge(int a, int b);
        void add_edge(int a, int b, double w);
        void set_weight(int a, int b, double w, bool new_edge = true);
        void transfer_weight(int a, int b, int c, int d, double w);
        bool allowed_weight(double w) const; //checks if a given weight is within the allowed bounds for this graph
        std::map<int,double> vertex_adjacencies(int v) const;
        const std::map<std::pair<int,int>, double>& edge_list_map() const;
        std::pair<std::pair<int,int>, std::pair<int,int> > sample_pair_edges();

};

Graph Graph_from_edge_list(Rcpp::NumericMatrix EdgeList, double lower_bound=0, double upper_bound=DBL_MAX, bool directed=false);


#endif /*GRAPH_H*/
