#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <cmath>
#include "ssmatrix.hpp"


void SSMatrix::empty_SSMatrix_initialization(Rcpp::IntegerMatrix M){
    n = M.nrow();
    m = M.ncol();
    B = std::vector<std::vector<bool> > (n, std::vector<bool>(m, 0));
    sampling_vector = SVector<std::pair<int,int>> ();
    nz = 0;
    nz_prob_numerator = 0;
    nz_rowsums = std::vector<int> (n, 0);
    nz_colsums = std::vector<int> (m, 0);
}



SSMatrix::SSMatrix(Rcpp::IntegerMatrix M, int min_row){
    empty_SSMatrix_initialization(M);
    nz_rowsums = std::vector<int> (n, 0);
    nz_colsums = std::vector<int> (m, 0);
    for (int i=min_row; i<n; ++i){
        for (int j=0; j<m; ++j){
            if (M(i,j) != 0){
                B[i][j] = true;
                sampling_vector.insert(std::make_pair(i, j));
                ++nz;
                ++nz_rowsums[i];
                ++nz_colsums[j];
            }
        }
    }
    nz_prob_numerator = pow(nz, 2) + nz;
    for (int i=min_row; i<n; ++i) 
        nz_prob_numerator -= pow(nz_rowsums[i], 2);
    for (int i=0; i<m; ++i) 
        nz_prob_numerator -= pow(nz_colsums[i], 2);
}

// SSMatrix::SSMatrix(Rcpp::IntegerMatrix M, Rcpp::LogicalMatrix sampling_area){
//     //unfinished
//     empty_SSMatrix_initialization(M);
//     for (int i=0; i<n; ++i){
//         for (int j=0; j<m; ++j){
//             if (sampling_area(i,j)){
//                 if (M(i,j)){
//                     sampling_vector.insert(std::make_pair(i, j));
//                     B[i][j]=true;
//                 }
//                 else{
//                     ++z;
//                     ++z_rowsums[i];
//                     ++z_colsums[j];
//                 }
//             }
//         }
//     }
//     //compute zero_prob_numerator here (!!!)
// }

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
        //Rcpp::Rcout << "a = (" << a.first << "," << a.second << "),  ";
        els[0] = a.first;
        els[2] = a.second;
        auto b = sample_element();
        //Rcpp::Rcout << "b = (" << b.first << "," << b.second << ")" << std::endl;
        els[1] = b.first;
        els[3] = b.second;
        if (els[0] != els[1] and els[2] != els[3])
            return els;
        ++k;
        if (k>1000*m*n){
            Rcpp::Rcout << "There might not be two elements satisfying the conditions" << std::endl;
            throw "There might not be two elements satisfying the conditions";
        }
    }
}

void SSMatrix::insert(std::pair<int,int> p){
    int x = p.first, y = p.second;
    if (not B[x][y]){
        sampling_vector.insert(p);
        B[x][y] = true;
        nz_prob_numerator += pow(nz_rowsums[x], 2) + pow(nz_colsums[y], 2);
        nz_prob_numerator += 2 * nz +2;
        ++nz_rowsums[x];
        ++nz_colsums[y];
        ++nz;
        nz_prob_numerator -= pow(nz_rowsums[x], 2) + pow(nz_colsums[y], 2);
    }
}

void SSMatrix::remove(std::pair<int,int> p){
    int x = p.first, y = p.second;
    if (B[x][y]){
        sampling_vector.remove(p);
        B[x][y] = false;
        nz_prob_numerator += pow(nz_rowsums[x], 2) + pow(nz_colsums[y], 2);
        nz_prob_numerator -= 2 * nz;
        --nz_rowsums[x];
        --nz_colsums[y];
        --nz;
        nz_prob_numerator -= pow(nz_rowsums[x], 2) + pow(nz_colsums[y], 2);
    }
}

int SSMatrix::sample_n_invariant_steps(){
    double p = nz_prob_numerator / (n * m * (n-1) * (m-1));
    // Rcpp::Rcout << "nz = " << nz << std::endl;
    // for (int i=0; i<n; ++i){
    //     Rcpp::Rcout << "nz_rowsums["<<i<<"] = " << nz_rowsums[i] << std::endl;
    // }
    // for (int i=0; i<m; ++i){
    //     Rcpp::Rcout << "nz_colsums["<<i<<"] = " << nz_colsums[i] << std::endl;
    // }
    // Rcpp::Rcout << "nz_prob_numerator = " << nz_prob_numerator << std::endl;
    // Rcpp::Rcout << "p = " << p << std::endl;
    Rcpp::NumericVector s = Rcpp::rgeom(1, p);
    return (int) s[0];
}





