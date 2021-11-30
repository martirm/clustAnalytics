#' clustAnalytics
#' 
#' This package evaluates the stability and significance of clusters in 
#' \code{igraph} graphs. Supports weighted and unweighted graphs.
#' 
#' Extensions to weighted graphs of multiple functions present in \code{igraph} 
#' are provided, such as scoring functions or edge rewiring methods.
#' 
#' @docType package
#' @encoding UTF-8
#' @author Martí Renedo Mirambell
#' @import Rcpp igraph
#' @importFrom Rcpp evalCpp
#' @importFrom Rdpack reprompt
#' @useDynLib clustAnalytics
#' @name clustAnalytics
NULL  