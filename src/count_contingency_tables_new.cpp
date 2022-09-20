#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include "ssmatrix.h"
#include "auxiliary_functions.h"

//using namespace Rcpp;


//' Performs a step of the Markov Chain Monte Carlo method
//' 
//' Modifies the matrix while keeping the column and row sums constant, as well as
//' leaving the positions strictly preceding (k,l) in lexicographical order invariant.
//' 
//' @param M matrix
//' @param k,l Coordinates of the first element that is not invariant 
//' @return boolean indicating whether the step left the matrix invariant
//' @keywords internal
// [[Rcpp::export]]  
bool walk_step(Rcpp::IntegerMatrix &M, int k, int l){
    //note that p1, p2, k, l are already 0-indexed
    int n=M.nrow(), m=M.ncol();
    
    //boundary conditions
    if (l>=m-1){
        return(walk_step(M,k+1,0));
    }
    if (k>=n-1){
        return false;
    }
    
    Rcpp::IntegerVector rows = Rcpp::sample(n - k, 2);
    Rcpp::IntegerVector cols = Rcpp::sample(m, 2);
    int i1 = rows[0] + k - 1, i2 = rows[1] + k - 1;
    int j1 = cols[0] - 1, j2 = cols[1] - 1;
    
    if (i1 == k or i2 == k){
        if (j1 < l or j2 < l){
            //these conditions would produce a change in the invariant part of the matrix.
            //In this case we simply resample i and j (equivalent to the recursive call)
            return walk_step(M, k, l);
        }
    }
    
    if(M(i1,j2) >= 1 and M(i2,j1)>=1){
        //std::cout << "(i1, j1)=(" << i1 << "," << j1 << ")" << std::endl;
        //std::cout << "(i2, j2)=(" << i2 << "," << j2 << ")" << std::endl;
        //print_Rcpp(M);
        --M(i1,j2);
        --M(i2,j1);
        ++M(i1,j1);
        ++M(i2,j2);
        return true;
    }
    return false;
}

void walk_n_steps(Rcpp::IntegerMatrix& M, int k, int l, int n){
    // repeats walk_step(M,k,l) until n steps are performed
    while (n > 0){
        walk_step(M, k, l); 
        --n;
    }
}

//' Estimates |H_i|/|H_{i+1}|
//' 
//' Estimates the fraction of elements of H_i that are also in H_{i+1} (where i=(k,l))
//' @param M matrix
//' @param k,l Coordinates of the first element that is not invariant 
//' @param error error for the convergence of the method
//' @return value of H_i/H_{i+1}
double estimate_H_fraction(const Rcpp::IntegerMatrix& M, int k, int l, double error=0.1){
    
    // boundary conditions (these fractions are 1 by definition of H)
    if (l==M.ncol()-1 or k==M.nrow()-1){
        return 1;
    }
    
    double H_i=0, H_i_plus=0; //we will determine the fraction of elements in H_i that are also in H_i_plus
    Rcpp::IntegerMatrix J(Rcpp::clone(M)); //this is necessary to avoid modifying M)
    int M_kl = M(k,l);
    //SSMatrix sampling_matrix(M, i);
    
    walk_n_steps(J, k, l, 1000);
    int iter = 0;
    double ratio = -1;
    bool is_H_i_plus = (J(k,l) == M_kl);
    while(1){ //replace this with convergence condition
        if (walk_step(J, k,l)){
            is_H_i_plus = (J(k,l) == M_kl);
        }
        if ( is_H_i_plus )
            ++H_i_plus;
        ++H_i;
        double new_ratio = H_i/H_i_plus;
        if (iter%100000 == 0){
            if (fabs(new_ratio - ratio) < error){ //replace with proper conditions depending on desired error
                //std::cout << iter << std::endl << H_i << std::endl << H_i_plus << std::endl;
                return new_ratio;
            }
            else
                ratio = new_ratio;
        }
        ++iter;
    }
    
}

//' Estimates |H_0|/|H_r*|
//' 
//' This is the total number of contingency tables (of the same margins as M) divided 
//' by the number that match M until the r-th row (included, 0-indexed). Note that
//' if r==0, this is always 1 by definition. 
//' @param M contingency table
//' @param r row index
//' @param error error for the convergence of the method
// [[Rcpp::export]]  
double estimate_H_fraction_r_rows(const Rcpp::IntegerMatrix& M, int r, double error=0.1){
    double fraction = 1;
    if (r > M.nrow()){
        //in this case, the fraction is the total number of contingency tables
        return(estimate_H_fraction_r_rows(M, M.nrow(), error));
    }
    for (int i=0; i<r; ++i){
        for (int j=0; j<M.ncol()-1; ++j){
            fraction *= estimate_H_fraction(M, i, j, error);
        }
    }
    return fraction;
}

//' Estimates |H_i|/|H_{i+1}| for the first r rows
//' 
//' The product of all these ratios is is the total number of contingency tables (of the same margins as M) divided 
//' by the number that match M until the r-th row (included, 0-indexed). 
//' @param M contingency table
//' @param r row index
//' @param error error for the convergence of the method
//' @return NumericVector containing all the ratios
//' @keywords internal
// [[Rcpp::export]]  
Rcpp::NumericVector estimate_H_fractions(const Rcpp::IntegerMatrix& M, int r, double error=0.1){
    if (r > M.nrow()){
        //in this case, the fraction is the total number of contingency tables
        r = M.nrow();
    }
    Rcpp::NumericVector H_fracs(r*(M.ncol()-1), 1.0);
    for (int i=0; i<r; ++i){
        for (int j=0; j<M.ncol()-1; ++j){
            H_fracs[i*(M.ncol()-1) + j] = estimate_H_fraction(M, i, j, error);
            //fraction *= estimate_H_fraction(M, i, j, error);
        }
    }
    return H_fracs;
}


//' Contingency table from membership vectors
//' 
//' Given two membership vectors, returns the corresponding contingency table.
//' we assume the labels are >=1 and numbered consecutively. If not consecutive
//' (some labels are unused) this implementation still works, but will be less
//' efficient.
//' @param c1,c2 membership vectors (integer values containing the index of each community)
//' @keywords internal
// [[Rcpp::export]]    
Rcpp::IntegerMatrix c_rs_table(const Rcpp::NumericVector& c1, const Rcpp::NumericVector& c2){
    
    int n = c1.size();
    int n_labels_c1 =  *std::max_element(c1.begin(), c1.end());
    int n_labels_c2 =  *std::max_element(c2.begin(), c2.end());
    Rcpp::IntegerMatrix c_rs(n_labels_c1, n_labels_c2); //filled with 0s by default
    for (int i=0; i<n; ++i){
        ++c_rs(c1[i]-1, c2[i]-1);
    }
    return c_rs;
}





