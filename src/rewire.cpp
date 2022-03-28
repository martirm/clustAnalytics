#include <Rcpp.h>
#include <string>
#include <map>
#include <algorithm>
#include "graph.h"

using namespace Rcpp;
using namespace std;


bool randomization_step(Graph &g, string weight_sel = "max_weight"){
    auto p = g.sample_pair_edges();
    int a = p.first.first, c = p.first.second,
        b = p.second.first, d = p.second.second;
    if (b==c or a==d or a==b or c==d){
        return false;
    }
    //Rcout << "a=" << a << ", b=" << b << ", c=" << c << ", d=" << d << endl;
    double wAC = g.get_weight(a, c),
           wAD = g.get_weight(a, d),
           wBC = g.get_weight(b, c),
           wBD = g.get_weight(b, d);

    double t;
    if (weight_sel == "max_weight"){
        double upper_bound = g.get_upper_bound();
        t = std::min({wAC, wBD, upper_bound-wAD, upper_bound-wBC});
        //Rcout << "t = " << t << endl;
    }
    else t = (wAC+wBD-wAD-wBC)/2;

    wAC -= t;
    wAD += t;
    wBC += t;
    wBD -= t;
    if ( g.allowed_weight(wAC) && g.allowed_weight(wAD) &&
         g.allowed_weight(wBC) && g.allowed_weight(wBD)){
        g.set_weight(a, c, wAC);
        g.set_weight(a, d, wAD);
        g.set_weight(b, c, wBC);
        g.set_weight(b, d, wBD);
        //Rcout << "switching" << endl;
        return true;
    }
    //Rcout << "not switching" << endl;
    return false;
}

void randomize_g(Graph &g, double Q, string weight_sel = "const_var"){
    int itermax = (int) (Q*g.get_size()), performed_switchings = 0;
    //Rcout << g.get_size() << endl;
    for (int i=0; i<itermax; ++i){
        //Rcout << i << endl;
        if (randomization_step(g, weight_sel)) {
            ++performed_switchings;
        }
    }
    // Rcout << "done" << endl;
    // Rcout << "Switchings performed: " << performed_switchings << "/"
    //       << itermax << endl;
}


//[[Rcpp::export]]
NumericMatrix randomize(NumericMatrix EdgeList, double Q, std::string weight_sel = "const_var",
                        double lower_bound=0, double upper_bound=-1, bool directed = false){
    if (upper_bound == -1) upper_bound = DBL_MAX;
    Graph g = Graph_from_edge_list(EdgeList, lower_bound, upper_bound, directed);
    randomize_g(g, Q, weight_sel);
    //Rcout << "converting back to edge list" << endl;
    return g.numericmatrix_edgelist();
}
