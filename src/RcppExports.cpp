// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// TPR_coms_Rcpp
NumericVector TPR_coms_Rcpp(const IntegerVector& triangles, const IntegerVector& com);
RcppExport SEXP _clustAnalytics_TPR_coms_Rcpp(SEXP trianglesSEXP, SEXP comSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type triangles(trianglesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type com(comSEXP);
    rcpp_result_gen = Rcpp::wrap(TPR_coms_Rcpp(triangles, com));
    return rcpp_result_gen;
END_RCPP
}
// triangle_participation_ratio_Rcpp
double triangle_participation_ratio_Rcpp(const NumericMatrix& EdgeList);
RcppExport SEXP _clustAnalytics_triangle_participation_ratio_Rcpp(SEXP EdgeListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type EdgeList(EdgeListSEXP);
    rcpp_result_gen = Rcpp::wrap(triangle_participation_ratio_Rcpp(EdgeList));
    return rcpp_result_gen;
END_RCPP
}
// walk_step
bool walk_step(IntegerMatrix& M, int min_row);
RcppExport SEXP _clustAnalytics_walk_step(SEXP MSEXP, SEXP min_rowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix& >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type min_row(min_rowSEXP);
    rcpp_result_gen = Rcpp::wrap(walk_step(M, min_row));
    return rcpp_result_gen;
END_RCPP
}
// walk_k_steps
void walk_k_steps(IntegerMatrix& M, int min_row, int k);
RcppExport SEXP _clustAnalytics_walk_k_steps(SEXP MSEXP, SEXP min_rowSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix& >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type min_row(min_rowSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    walk_k_steps(M, min_row, k);
    return R_NilValue;
END_RCPP
}
// sample_fraction_H_i
int sample_fraction_H_i(IntegerMatrix M, int i, double error);
RcppExport SEXP _clustAnalytics_sample_fraction_H_i(SEXP MSEXP, SEXP iSEXP, SEXP errorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< double >::type error(errorSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_fraction_H_i(M, i, error));
    return rcpp_result_gen;
END_RCPP
}
// c_rs_table
IntegerMatrix c_rs_table(const NumericVector& c1, const NumericVector& c2);
RcppExport SEXP _clustAnalytics_c_rs_table(SEXP c1SEXP, SEXP c2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type c1(c1SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type c2(c2SEXP);
    rcpp_result_gen = Rcpp::wrap(c_rs_table(c1, c2));
    return rcpp_result_gen;
END_RCPP
}
// count_contingency_tables
int count_contingency_tables(const NumericVector& c1, const NumericVector& c2, double error);
RcppExport SEXP _clustAnalytics_count_contingency_tables(SEXP c1SEXP, SEXP c2SEXP, SEXP errorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type c1(c1SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type c2(c2SEXP);
    Rcpp::traits::input_parameter< double >::type error(errorSEXP);
    rcpp_result_gen = Rcpp::wrap(count_contingency_tables(c1, c2, error));
    return rcpp_result_gen;
END_RCPP
}
// count_labels
IntegerVector count_labels(const IntegerVector& c);
RcppExport SEXP _clustAnalytics_count_labels(SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(count_labels(c));
    return rcpp_result_gen;
END_RCPP
}
// mutual_information_Cpp
double mutual_information_Cpp(const IntegerVector& c1, const IntegerVector& c2, const IntegerVector& a, const IntegerVector& b);
RcppExport SEXP _clustAnalytics_mutual_information_Cpp(SEXP c1SEXP, SEXP c2SEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type c1(c1SEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type c2(c2SEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(mutual_information_Cpp(c1, c2, a, b));
    return rcpp_result_gen;
END_RCPP
}
// vector_c_rs
IntegerVector vector_c_rs(const IntegerVector& c1, const IntegerVector& c2);
RcppExport SEXP _clustAnalytics_vector_c_rs(SEXP c1SEXP, SEXP c2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type c1(c1SEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type c2(c2SEXP);
    rcpp_result_gen = Rcpp::wrap(vector_c_rs(c1, c2));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _clustAnalytics_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// resampled_edgelist
NumericMatrix resampled_edgelist(NumericMatrix el, NumericVector s);
RcppExport SEXP _clustAnalytics_resampled_edgelist(SEXP elSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type el(elSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(resampled_edgelist(el, s));
    return rcpp_result_gen;
END_RCPP
}
// randomize
NumericMatrix randomize(NumericMatrix EdgeList, double Q, std::string weight_sel, double lower_bound, double upper_bound);
RcppExport SEXP _clustAnalytics_randomize(SEXP EdgeListSEXP, SEXP QSEXP, SEXP weight_selSEXP, SEXP lower_boundSEXP, SEXP upper_boundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type EdgeList(EdgeListSEXP);
    Rcpp::traits::input_parameter< double >::type Q(QSEXP);
    Rcpp::traits::input_parameter< std::string >::type weight_sel(weight_selSEXP);
    Rcpp::traits::input_parameter< double >::type lower_bound(lower_boundSEXP);
    Rcpp::traits::input_parameter< double >::type upper_bound(upper_boundSEXP);
    rcpp_result_gen = Rcpp::wrap(randomize(EdgeList, Q, weight_sel, lower_bound, upper_bound));
    return rcpp_result_gen;
END_RCPP
}
// cluster_auxiliary_values_Rcpp
NumericMatrix cluster_auxiliary_values_Rcpp(const NumericMatrix& EdgeList, const IntegerVector& memb);
RcppExport SEXP _clustAnalytics_cluster_auxiliary_values_Rcpp(SEXP EdgeListSEXP, SEXP membSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type EdgeList(EdgeListSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type memb(membSEXP);
    rcpp_result_gen = Rcpp::wrap(cluster_auxiliary_values_Rcpp(EdgeList, memb));
    return rcpp_result_gen;
END_RCPP
}
// density_ratio_Rcpp
double density_ratio_Rcpp(const NumericMatrix& aux_vals);
RcppExport SEXP _clustAnalytics_density_ratio_Rcpp(SEXP aux_valsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type aux_vals(aux_valsSEXP);
    rcpp_result_gen = Rcpp::wrap(density_ratio_Rcpp(aux_vals));
    return rcpp_result_gen;
END_RCPP
}
// local_density_ratio_Rcpp
NumericVector local_density_ratio_Rcpp(const NumericMatrix& aux_vals);
RcppExport SEXP _clustAnalytics_local_density_ratio_Rcpp(SEXP aux_valsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type aux_vals(aux_valsSEXP);
    rcpp_result_gen = Rcpp::wrap(local_density_ratio_Rcpp(aux_vals));
    return rcpp_result_gen;
END_RCPP
}
// FOMD_Rcpp
NumericVector FOMD_Rcpp(const NumericMatrix& EdgeList, const IntegerVector& memb);
RcppExport SEXP _clustAnalytics_FOMD_Rcpp(SEXP EdgeListSEXP, SEXP membSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type EdgeList(EdgeListSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type memb(membSEXP);
    rcpp_result_gen = Rcpp::wrap(FOMD_Rcpp(EdgeList, memb));
    return rcpp_result_gen;
END_RCPP
}
// out_degree_fractions_Rcpp
NumericMatrix out_degree_fractions_Rcpp(const NumericMatrix& EdgeList, const IntegerVector& memb);
RcppExport SEXP _clustAnalytics_out_degree_fractions_Rcpp(SEXP EdgeListSEXP, SEXP membSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type EdgeList(EdgeListSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type memb(membSEXP);
    rcpp_result_gen = Rcpp::wrap(out_degree_fractions_Rcpp(EdgeList, memb));
    return rcpp_result_gen;
END_RCPP
}
// clustering_coefficient_Rcpp
double clustering_coefficient_Rcpp(const NumericMatrix& EdgeList, double lower_bound, double upper_bound);
RcppExport SEXP _clustAnalytics_clustering_coefficient_Rcpp(SEXP EdgeListSEXP, SEXP lower_boundSEXP, SEXP upper_boundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type EdgeList(EdgeListSEXP);
    Rcpp::traits::input_parameter< double >::type lower_bound(lower_boundSEXP);
    Rcpp::traits::input_parameter< double >::type upper_bound(upper_boundSEXP);
    rcpp_result_gen = Rcpp::wrap(clustering_coefficient_Rcpp(EdgeList, lower_bound, upper_bound));
    return rcpp_result_gen;
END_RCPP
}
// transitivity_Rcpp
double transitivity_Rcpp(const NumericMatrix& EdgeList, double lower_bound, double upper_bound);
RcppExport SEXP _clustAnalytics_transitivity_Rcpp(SEXP EdgeListSEXP, SEXP lower_boundSEXP, SEXP upper_boundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type EdgeList(EdgeListSEXP);
    Rcpp::traits::input_parameter< double >::type lower_bound(lower_boundSEXP);
    Rcpp::traits::input_parameter< double >::type upper_bound(upper_boundSEXP);
    rcpp_result_gen = Rcpp::wrap(transitivity_Rcpp(EdgeList, lower_bound, upper_bound));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_clustAnalytics_TPR_coms_Rcpp", (DL_FUNC) &_clustAnalytics_TPR_coms_Rcpp, 2},
    {"_clustAnalytics_triangle_participation_ratio_Rcpp", (DL_FUNC) &_clustAnalytics_triangle_participation_ratio_Rcpp, 1},
    {"_clustAnalytics_walk_step", (DL_FUNC) &_clustAnalytics_walk_step, 2},
    {"_clustAnalytics_walk_k_steps", (DL_FUNC) &_clustAnalytics_walk_k_steps, 3},
    {"_clustAnalytics_sample_fraction_H_i", (DL_FUNC) &_clustAnalytics_sample_fraction_H_i, 3},
    {"_clustAnalytics_c_rs_table", (DL_FUNC) &_clustAnalytics_c_rs_table, 2},
    {"_clustAnalytics_count_contingency_tables", (DL_FUNC) &_clustAnalytics_count_contingency_tables, 3},
    {"_clustAnalytics_count_labels", (DL_FUNC) &_clustAnalytics_count_labels, 1},
    {"_clustAnalytics_mutual_information_Cpp", (DL_FUNC) &_clustAnalytics_mutual_information_Cpp, 4},
    {"_clustAnalytics_vector_c_rs", (DL_FUNC) &_clustAnalytics_vector_c_rs, 2},
    {"_clustAnalytics_rcpp_hello_world", (DL_FUNC) &_clustAnalytics_rcpp_hello_world, 0},
    {"_clustAnalytics_resampled_edgelist", (DL_FUNC) &_clustAnalytics_resampled_edgelist, 2},
    {"_clustAnalytics_randomize", (DL_FUNC) &_clustAnalytics_randomize, 5},
    {"_clustAnalytics_cluster_auxiliary_values_Rcpp", (DL_FUNC) &_clustAnalytics_cluster_auxiliary_values_Rcpp, 2},
    {"_clustAnalytics_density_ratio_Rcpp", (DL_FUNC) &_clustAnalytics_density_ratio_Rcpp, 1},
    {"_clustAnalytics_local_density_ratio_Rcpp", (DL_FUNC) &_clustAnalytics_local_density_ratio_Rcpp, 1},
    {"_clustAnalytics_FOMD_Rcpp", (DL_FUNC) &_clustAnalytics_FOMD_Rcpp, 2},
    {"_clustAnalytics_out_degree_fractions_Rcpp", (DL_FUNC) &_clustAnalytics_out_degree_fractions_Rcpp, 2},
    {"_clustAnalytics_clustering_coefficient_Rcpp", (DL_FUNC) &_clustAnalytics_clustering_coefficient_Rcpp, 3},
    {"_clustAnalytics_transitivity_Rcpp", (DL_FUNC) &_clustAnalytics_transitivity_Rcpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_clustAnalytics(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
