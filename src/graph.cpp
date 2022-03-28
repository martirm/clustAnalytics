#include <vector>
#include <map>
#include <unordered_map>
#include <cfloat>
#include "graph.h"
#include "cantor_hash.h"

using namespace std;


Graph::Graph( vector<Edge> v, int graph_order, double lb, double ub, bool d){
    order = graph_order;
    lower_bound = lb;
    upper_bound = ub;
    directed = d;
    edge_list_m = unordered_map<pair<int,int>, double, CantorHash>();
    adjacencies_list = vector<map<int, double> > (order);
    size=v.size();
    vector<pair<int,int> > edge_vector(size); //same as v but without weight information. Needed to call the sampling vector constructor.
    for (int i=0; i<size; ++i){
        pair<int,int> ab = make_pair(v[i].a, v[i].b);
        edge_list_m[ab] = v[i].weight;
        edge_vector[i] = ab;
        adjacencies_list[v[i].a] [v[i].b] = v[i].weight;
        adjacencies_list[v[i].b] [v[i].a] = v[i].weight;
    }
    sampling_vector = SVector<pair<int,int> >(edge_vector);
}

//wrapper for creating Graph object from a NumericMatrix edgelist
Graph Graph_from_edge_list(Rcpp::NumericMatrix EdgeList, double lower_bound,
                           double upper_bound, bool directed){
    int max_index = 0;
    int n = EdgeList.nrow();
    vector<Edge> v(n);
    for (int i=0; i<n; ++i){
        //substract 1 to convert to 0-based indices
        v[i] = {int(EdgeList(i,0))-1, int(EdgeList(i,1))-1, EdgeList(i,2)};
        if (v[i].b > max_index) max_index = v[i].b;
    }
    return Graph(v, max_index+1, lower_bound, upper_bound, directed);
}

int Graph::get_size() const{
    return size;
}

int Graph::get_order() const{
    return order;
}

double Graph::get_upper_bound() const{
    return upper_bound;
}

double Graph::get_lower_bound() const{
    return lower_bound;
}

void Graph::update_size(){
    size = edge_list_m.size();
}

double Graph::get_weight(int a, int b) const{
    if (directed and b<a) swap(a, b);
    if (b>=order) {
        return -1;
    }
    pair<int,int> ab = make_pair(a, b);
    auto it = edge_list_m.find(ab);
    if (it == edge_list_m.end()) return 0;
    else return it->second;
}

bool Graph::adjacent(int a, int b) const{
    if (directed and b<a) swap(a, b);
    return edge_list_m.find(make_pair(a, b)) == edge_list_m.end();
}

void Graph::delete_edge(int a, int b){
    if (directed and b<a) swap(a, b);
    pair<int,int> ab = make_pair(a, b);
    edge_list_m.erase(ab);
    adjacencies_list[a].erase(b);
    adjacencies_list[b].erase(a);
    sampling_vector.remove(ab);
    update_size();
}

void Graph::add_edge(int a, int b, double w){
    set_weight(a, b, w, true);
    update_size();
}


void Graph::set_weight(int a, int b, double w, bool new_edge){
    if (w==0){
        delete_edge(a, b);
        return;
    }
    if (directed and b<a) swap(a, b);
    pair<int,int> ab = make_pair(a, b);
    edge_list_m[ab] = w;
    adjacencies_list[a][b] = w;
    adjacencies_list[b][a] = w;
    if (new_edge) sampling_vector.insert(ab);
    update_size();
}


//samples pair of edges
pair<pair<int,int>, pair<int,int> > Graph::sample_pair_edges(){
    pair<int,int> ac = sampling_vector.rand_el(), bd;
    bool same = true;
    while (same){
        bd = sampling_vector.rand_el();
        if (bd != ac){
            same = false;
        }
    }
    return make_pair(ac, bd);
}

bool Graph::allowed_weight(double w) const{
    return (w <= upper_bound && w >= lower_bound);
}

std::map<int,double> Graph::vertex_adjacencies(int v) const{
    return adjacencies_list[v];
}

const std::unordered_map<std::pair<int,int>, double, CantorHash>& Graph::edge_list_map() const{
    return edge_list_m;
}

//returns edge list as a NumericMatrix RCpp object
Rcpp::NumericMatrix Graph::numericmatrix_edgelist(){
    Rcpp::NumericMatrix el_nm(size, 3);
    int i=0;
    for (auto it=edge_list_m.begin(); it!=edge_list_m.end(); ++it){
        //sum 1 to convert back to 1-based indices
        el_nm(i,0) = it->first.first+1;
        el_nm(i,1) = it->first.second+1;
        el_nm(i,2) = it->second;
        ++i;
    }
    return el_nm;
}
