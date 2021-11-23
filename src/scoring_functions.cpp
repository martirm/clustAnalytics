#include <Rcpp.h>
#include <map>
#include <vector>
#include "graph.h"
using namespace Rcpp;

struct cluster_values{
    double m;
    double c;
    int n;
};

//maybe this one won't be necessary
std::vector<cluster_values> cluster_auxiliary_values(const Graph& g, const std::vector<int>& memb){
    // memb contains the indices of the communities. We assume they are indexed from 1 to n_com
    int n_com = *std::max_element(memb.begin(), memb.end());
    std::vector<cluster_values> v(n_com, {0,0,0});
    auto edge_list_map = g.edge_list_map();
    for (auto it=edge_list_map.begin(); it!=edge_list_map.end(); ++it){
        int a = it->first.first, b = it->first.second;

        if (memb[a] == memb[b]){
            //intra community edge
            v[memb[a]].m += it->second;
        }
        else{
            //inter community edge
            v[memb[a]].c += it->second;
            v[memb[b]].c += it->second;
        }
    }
    for (unsigned int i = 0; i < memb.size(); ++i){
        v[memb[i]].n++;
    }

    return v;
}

//[[Rcpp::export]]
NumericMatrix cluster_auxiliary_values_Rcpp(const NumericMatrix& EdgeList, const IntegerVector& memb){
    /* memb contains the indices of the communities. We assume they are indexed from 1 to n_com
    returns a numeric matrix where rows are communities, and cols are,
    respectively, m_w, n, c_w. */
    int n_com = *std::max_element(memb.begin(), memb.end());
    NumericMatrix values(n_com, 3); //filled with 0s by default
    for (int i = 0; i < EdgeList.nrow(); ++i){
        // we convert the vertex indices to 0-based
        int a = EdgeList(i, 0)-1, b = EdgeList(i, 1)-1;

        if (memb[a] == memb[b]){
            //intra community edge
            /* indices of communities in memb are 1-based, but the Numericmatrix
            has 0-based indices in C++, so we substract 1*/
            values(memb[a]-1, 0) += EdgeList(i, 2);
        }
        else{
            //inter community edge
            values(memb[a]-1, 2) += EdgeList(i, 2);
            values(memb[b]-1, 2) += EdgeList(i, 2);
        }
    }
    for (int i = 0; i < memb.size(); ++i){
        values(memb[i]-1,1)++;
    }

    return values;
}





//[[Rcpp::export]]
double density_ratio_Rcpp(const NumericMatrix& aux_vals){
    double internal_weight=0, external_weight=0;
    int potential_internal=0, potential_external=0; //potential number of edges of each kind
    int n=0, n_row = aux_vals.nrow();
    for (int i=0; i<n_row; ++i) n += aux_vals(i, 1);  //total number of vertices
    for (int i=0; i<n_row; ++i){
        internal_weight += aux_vals(i, 0);
        external_weight += aux_vals(i, 2);
        potential_internal += aux_vals(i, 1) * (aux_vals(i, 1) -1);
        potential_external += aux_vals(i, 1) * (n - aux_vals(i, 1));
    }
    //remove duplicates:
    potential_internal /= 2;
    potential_external /= 2;

    return 1. - ((external_weight/potential_external) / (internal_weight/potential_internal));
}

IntegerVector potential_external_edges(const NumericVector& com_sizes){
    int n=0, n_coms = com_sizes.size();
    IntegerVector pe (n_coms, 0);
    for (int i=0; i<n_coms; ++i) n += com_sizes[i];
    for (int i=0; i<n_coms; ++i){
        pe[i] += com_sizes[i] * (n - com_sizes[i]);
    }
    return pe;
}

//[[Rcpp::export]]
NumericVector local_density_ratio_Rcpp(const NumericMatrix& aux_vals){
    NumericVector pi = aux_vals(_, 1) * (aux_vals(_, 1) - 1) / 2;
    NumericVector pe = as<NumericVector>(potential_external_edges(aux_vals(_, 1)));
    return 1. - ((aux_vals(_,2) / pe) / (aux_vals(_, 0) / pi));
}

//[[Rcpp::export]]
NumericVector FOMD_Rcpp(const NumericMatrix& EdgeList, const IntegerVector& memb){
    // Inputs are the same as above. Returns a vector with the FOMD (fraction
    // over median degree)
    int n_com = *std::max_element(memb.begin(), memb.end());
    int n = memb.size();
    std::vector<double> degree(n, 0), internal_degree(n, 0);
    for (int i = 0; i < EdgeList.nrow(); ++i){
        // correcting indices
        int a = EdgeList(i, 0)-1, b = EdgeList(i, 1)-1;

        if (memb[a] == memb[b]){
            internal_degree[a] += EdgeList(i, 2);
            internal_degree[b] += EdgeList(i, 2);
        }
        degree[a] += EdgeList(i, 2);
        degree[b] += EdgeList(i, 2);
    }

    // computing the median degree
    std::nth_element(degree.begin(), degree.begin() + n/2, degree.end());
    double median = degree[n/2];
    if (n%2 != 0){
        std::nth_element(degree.begin(), degree.begin() + n/2 - 1, degree.end());
        median = (median + degree[n/2 - 1]) / 2.0;
    }

    NumericVector FOMD(n_com, 0.0), com_size(n, 0.0);
    for (int i=0; i<n; ++i){
        int com = memb[i] - 1; // 0-based index
        ++com_size[com];
        if (internal_degree[i] > median) ++FOMD[com];
    }
    FOMD = FOMD / com_size;

    return FOMD;
}

//[[Rcpp::export]]
NumericMatrix out_degree_fractions_Rcpp(const NumericMatrix& EdgeList, const IntegerVector& memb){
    /* Inputs are the same as above. Returns a matrix where rows are communities
    and cols are maximum-ODF, average-ODF and flake-ODF respectively.
    We assume the graph is undirected (directed version can be easily implemented
    if needed)*/
    int n_com = *std::max_element(memb.begin(), memb.end());
    int n = memb.size();
    std::vector<double> degree(n, 0), out_degree(n, 0), com_size(n_com, 0);
    NumericMatrix values(n_com, 3);

    //read edge list to compute degree and out degree of all vertices.
    for (int i = 0; i < EdgeList.nrow(); ++i){
        //correcting indices (same as above)
        int a = EdgeList(i, 0)-1, b = EdgeList(i, 1)-1;

        if (memb[a] != memb[b]){
            out_degree[a] += EdgeList(i, 2);
            out_degree[b] += EdgeList(i, 2);
        }
        degree[a] += EdgeList(i, 2);
        degree[b] += EdgeList(i, 2);
    }

    //now use them to compute each type of ODF
    for(int i=0; i<n; ++i){
        double odf = out_degree[i] / degree[i];
        int com = memb[i] - 1;
        //max-ODF
        if (odf > values(com, 0)){
            values(com,0) = odf; //update  if necessary
        }
        //average-ODF (for now we just sum all odf's to compute the mean later)
        values(com, 1) += odf;
        com_size[com] += 1;
        //flake-ODF
        if (odf > 0.5){
            values(com, 2)++;
        }
    }
    //we divide the sum of odf's of each community by its size to get the mean
    //same for flake-ODF to get the fraction of nodes
    for (int i=0; i<n_com; ++i){
        values(i, 1) /= com_size[i];
        values(i, 2) /= com_size[i];
    }
    return values;
}
