#include <Rcpp.h>
#include <set>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix resampled_edgelist(NumericMatrix el, NumericVector s) {
    //note that the edge indices are converted to 0 based for calculations here and then returned to 1 based
    std::vector<std::set<int> > inv_s(s.size());
    for (int i=0; i<(int)s.size(); ++i) inv_s[s[i]-1].insert(i);
    std::vector<std::pair<int,int> > el_rs; //resampled edge list
    std::vector<double> weights_rs; //resampled edge weights
    for (int i=0; i<el.nrow(); ++i){
        int e1=el(i,0)-1, e2=el(i,1)-1;
        for (auto j = inv_s[e1].begin(); j != inv_s[e1].end(); ++j){
            for (auto k = inv_s[e2].begin(); k != inv_s[e2].end(); ++k){
                el_rs.push_back(std::make_pair(*j,*k));
                weights_rs.push_back(el(i,2));
            }
        }
    }
    NumericMatrix weighted_el_rs (weights_rs.size(), 3);
    for (int i=0; i<(int)weights_rs.size(); ++i){ 
        weighted_el_rs(i,0) = el_rs[i].first+1;
        weighted_el_rs(i,1) = el_rs[i].second+1;
        weighted_el_rs(i,2) = weights_rs[i];
    }
    return weighted_el_rs;
    
}