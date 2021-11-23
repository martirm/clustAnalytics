#ifndef SSMATRIX_HPP // include guard
#define SSMATRIX_HPP

#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include "svector.h"

// SSMatrix stands for Sparse Sampling Matrix
class SSMatrix{
    // contains information on which elements of a matrix are non zero, and allows
    // uniform sampling. Only samples the indices, not the elements themselves.
        SVector<std::pair<int,int> > sampling_vector;
        int n, m, z;  //z is the number of zeroes
        std::vector<int> z_rowsums, z_colsums;
        double zero_prob_numerator;
        std::vector<std::vector<bool>> B; //used to keep track of the zeroes;
    public:
        void empty_SSMatrix_initialization(Rcpp::IntegerMatrix M);
        SSMatrix(Rcpp::IntegerMatrix M, int min_row);
        SSMatrix(Rcpp::IntegerMatrix M) : SSMatrix(M, 0) {};
        SSMatrix(Rcpp::IntegerMatrix M, Rcpp::LogicalMatrix sampling_area);
        
        std::pair<int, int> sample_element();
        std::vector<int> sample_rw_step();
        
        void insert(std::pair<int,int> p);
        void insert(int x, int y) {insert(std::make_pair(x, y));};

        void remove(std::pair<int,int> p);
        void remove(int x, int y) {remove(std::make_pair(x, y));};
        
        int sample_n_invariant_steps();
        
        
};


#endif /* SSMATRIX_HPP */

