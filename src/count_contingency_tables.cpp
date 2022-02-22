#include <Rcpp.h>
#include <cmath>
#include "ssmatrix.h"

using namespace Rcpp;


// [[Rcpp::export]]  
bool walk_step(IntegerMatrix& M, int min_row){
    // samples rows and columns and performs the step.
    // returns true if successful, returns false if sampled step can't be performed
    //min_col is 0 indexed
    
    //sample i1,i2 from min_row:(nrow-1) and j1,j2 from 0:(ncol-1)
    IntegerVector rows = Rcpp::sample(M.nrow() - min_row, 2);
    IntegerVector cols = Rcpp::sample(M.ncol(), 2);
    int i1 = rows[0] + min_row - 1, i2 = rows[1] + min_row - 1;
    int j1 = cols[0] - 1, j2 = cols[1] - 1;
    
    if(M(i1,j2) >= 1 and M(i2,j1)>=1){
        --M(i1,j2);
        --M(i2,j1);
        ++M(i1,j1);
        ++M(i2,j2);
        return true;
    }
    return false;
}

void walk_effective_step(IntegerMatrix& M, SSMatrix& SM){
    std::vector<int> v = SM.sample_rw_step();
    int i1 = v[0], i2 = v[1], j1 = v[2], j2 = v[3];
    if (--M(i1, j1) == 0) SM.remove(i1, j1);
    if (--M(i2, j2) == 0) SM.remove(i2, j2);
    if (++M(i1, j2) == 1) SM.insert(i1, j2);
    if (++M(i2, j1) == 1) SM.insert(i2, j1);
}

// [[Rcpp::export]]  
void walk_k_steps(IntegerMatrix& M, int min_row, int k){
    // repeats walk_step until k steps are performed
    while (k > 0){
        walk_step(M, min_row); 
        --k;
    }
}

bool row_is_equal (const IntegerMatrix& M, const IntegerVector& v, int row){
    int n=v.size();
    for (int i=0; i<n; ++i){
        if (M(row,i) != v[i])
            return false;
    }
    return true;
}

// [[Rcpp::export]]  
int sample_fraction_H_i(IntegerMatrix M, int i, double error=0.1){
    double H_i=0, H_i_plus=0; //we will determine the fraction of elements in H_i that are also in H_i_plus
    IntegerMatrix J(Rcpp::clone(M)); //this is necessary to avoid modifying M)
    IntegerVector row_i = M(i,_);
    SSMatrix sampling_matrix(M, i);
    
    
    walk_k_steps(J, i, 1000);
    int iter=0;
    double ratio = -1;
    bool is_H_i_plus = false;
    while(1){ //replace this with convergence condition
        if (walk_step(J, i)){
            is_H_i_plus = row_is_equal(J, row_i, i);
        }
        if ( is_H_i_plus )
            ++H_i_plus;
        ++H_i;
        double new_ratio = H_i/H_i_plus;
        if (iter%100000 == 0){
            if (fabs(new_ratio - ratio) < error){ //replace with proper conditions depending on desired error
                //std::cout << iter << std::endl << H_i << std::endl << H_i_plus << std::endl;
                return round(new_ratio);
            }
            else
                ratio = new_ratio;
        }
        ++iter;
    }
}




// [[Rcpp::export]]    
IntegerMatrix c_rs_table(const NumericVector& c1, const NumericVector& c2){
    // we assume the labels are >=1 and numbered consecutively. If not consecutive
    // (some labels are unused) this implementation still works, but will be less
    // efficient.
    
    int n = c1.size();
    int n_labels_c1 =  *std::max_element(c1.begin(), c1.end());
    int n_labels_c2 =  *std::max_element(c2.begin(), c2.end());
    IntegerMatrix c_rs(n_labels_c1, n_labels_c2); //filled with 0s by default
    for (int i=0; i<n; ++i){
        ++c_rs(c1[i]-1, c2[i]-1);
    }
    return c_rs;
}
    
// [[Rcpp::export]]
int count_contingency_tables(const NumericVector& c1, const NumericVector& c2, double error=0.1){
    
    
    IntegerMatrix M = c_rs_table(c1,c2);
    int nrow = M.nrow(), n_contingency_tables = 1;
    if (nrow==1 or M.ncol()==1){
        return 1;
    }
    
    IntegerVector v_H(nrow-1);
    for (int i=0; i<nrow-1; ++i){
        v_H[i] = sample_fraction_H_i(M, i, error);
        n_contingency_tables *= v_H[i];
    }
    
    return n_contingency_tables;
}

