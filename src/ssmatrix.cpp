#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <cmath>
#include "ssmatrix.h"


void SSMatrix::empty_SSMatrix_initialization(Rcpp::IntegerMatrix M){
    n = M.nrow();
    m = M.ncol();
    B = std::vector<std::vector<bool> > (n, std::vector<bool>(m, 0));
    sampling_vector = SVector<std::pair<int,int>> ();
    z = 0;
    z_rowsums = std::vector<int> (n, 0);
    z_colsums = std::vector<int> (m, 0);
}

SSMatrix::SSMatrix(Rcpp::IntegerMatrix M, int min_row){
    empty_SSMatrix_initialization(M);
    z = (n-min_row) * m;
    z_rowsums = std::vector<int> (n, m);
    z_colsums = std::vector<int> (m, n);
    for (int i=min_row; i<n; ++i){
        for (int j=0; j<m; ++j){
            if (M(i,j) != 0){
                sampling_vector.insert(std::make_pair(i, j));
                B[i][j]=true;
                --z;
                --z_rowsums[i];
                --z_colsums[j];
            }
        }
    }
}

SSMatrix::SSMatrix(Rcpp::IntegerMatrix M, Rcpp::LogicalMatrix sampling_area){
    empty_SSMatrix_initialization(M);
    for (int i=0; i<n; ++i){
        for (int j=0; j<m; ++j){
            if (sampling_area(i,j)){
                if (M(i,j)){
                    sampling_vector.insert(std::make_pair(i, j));
                    B[i][j]=true;
                }
                else{
                    ++z;
                    ++z_rowsums[i];
                    ++z_colsums[j];
                }
            }
        }
    }
}

std::pair<int, int> SSMatrix::sample_element(){
    return sampling_vector.rand_el();
}

std::vector<int> SSMatrix::sample_rw_step(){
    // samples two elements that serve as a step for the random walk for the contingency table monte carlo sampling
    // That is, two elements that don't share neither row nor column
    // returned as the vector {i, i', j, j'}
    std::vector<int> els(4);
    int k=0;
    while (1){
        auto a = sample_element();
        els[0] = a.first;
        els[2] = a.second;
        auto b = sample_element();
        els[1] = b.first;
        els[3] = b.second;
        if (els[0] != els[1] and els[2] != els[3])
            return els;
        ++k;
        if (k>10*m*n){
            throw "There might not be two elements satisfying the conditions";
        }
    }
}

void SSMatrix::insert(std::pair<int,int> p){
    int x = p.first, y = p.second;
    if (not B[x][y]){
        sampling_vector.insert(p);
        B[x][y] = true;
        zero_prob_numerator += pow(z_rowsums[x], 2) + pow(z_colsums[y], 2);
        zero_prob_numerator -= 2 * z - 1;
        --z_rowsums[x];
        --z_colsums[y];
        --z;
        zero_prob_numerator -= pow(z_rowsums[x], 2) + pow(z_colsums[y], 2);
    }
}

void SSMatrix::remove(std::pair<int,int> p){
    int x = p.first, y = p.second;
    if (B[x][y]){
        sampling_vector.remove(p);
        B[x][y] = true;
        zero_prob_numerator += pow(z_rowsums[x], 2) + pow(z_colsums[y], 2);
        zero_prob_numerator += 2 * z + 1;
        ++z_rowsums[x];
        ++z_colsums[y];
        ++z;
        zero_prob_numerator -= pow(z_rowsums[x], 2) + pow(z_colsums[y], 2);
    }
}

int SSMatrix::sample_n_invariant_steps(){
    double p = zero_prob_numerator / (n * m * (n-1) * (m-1));
    Rcpp::NumericVector s = Rcpp::rgeom(1, p);
    return (int) s[0];
}





