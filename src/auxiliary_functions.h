#ifndef AUXILIARY_FUNCTIONS_H // include guard
#define AUXILIARY_FUNCTIONS_H

#include <Rcpp.h>

// [[Rcpp::export]]  
void print_Rcpp(Rcpp::IntegerMatrix M){
    int n=M.nrow(), m=M.ncol();
    for (int i=0; i<n; ++i){
        Rcpp::Rcout << (i==0 ? "(" : " ");
        for (int j=0; j<m; ++j){
            Rcpp::Rcout << " " << M(i,j) ;
        }
        if (i==n-1) Rcpp::Rcout << " )";
        Rcpp::Rcout << std::endl;
    }
}


#endif /*UXILIARY_FUNCTIONS_H*/