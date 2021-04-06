#include <Rcpp.h>
#include <vector>
#include <map>

using namespace Rcpp;

// IntegerVector is used to avoid copying vectors from R, but used just like an std::vector<int>

//[[Rcpp::export]]
IntegerVector count_labels (const IntegerVector& c){
    // Returns a vector of frequencies of vector c
    // VERY IMPORTANT: the returned vector is 1-indexed. That is, the first
    //      position (0 on the Rcpp side but 1 on the R side) contains the
    //      number of ones in c.

    int n = c.size();
    int m = *std::max_element(c.begin(), c.end());
    IntegerVector v(m, 0);
    for (int i=0; i<n; ++i) v[c[i]-1]++;
    return v;
}

std::map<std::pair<int,int>, int> count_c_rs(const IntegerVector& c1, 
                                             const IntegerVector& c2 ){

    int n = c1.size();
    std::map<std::pair<int,int>, int> c_rs;
    for (int i=0; i<n; ++i){
        std::pair<int,int> p (c1[i], c2[i]);
        auto em = c_rs.emplace(std::make_pair(p, 1));
        if (em.second == false){
            em.first->second++;
        }
    }
    return c_rs;
}

//[[Rcpp::export]]
double mutual_information_Cpp(const IntegerVector& c1, const IntegerVector& c2,
                               const IntegerVector& a, const IntegerVector& b){

    auto c_rs = count_c_rs(c1, c2);
    double n = c1.size();
    double MI = 0;
    for (auto const& x: c_rs){
        const int &r = x.first.first, s = x.first.second;
        double c = x.second, P_rs = c / n;
        //Rcout << "r: " << r << ",  s: " << s << std::endl;
        //Rcout << "P_rs: " << P_rs << ",  P_r: " << a[r-1]/n << ",  P_s: " << b[s-1]/n << std::endl;
        //Rcout << P_rs * n * n / double(a[r-1]) / double(b[s-1]) << std::endl << "************" << std::endl;
        MI += P_rs * log(P_rs * n * n / double(a[r-1]) / double(b[s-1]));
    }
    return MI;
}

//[[Rcpp::export]]
IntegerVector vector_c_rs(const IntegerVector& c1, const IntegerVector& c2){
    //returns a vector with al c_rs values
    auto c_rs = count_c_rs(c1, c2);
    int l = c_rs.size();
    IntegerVector v(l);
    int i=0;
    for (auto const& x: c_rs){
        v[i++] = x.second;
    }
    return v;
}

