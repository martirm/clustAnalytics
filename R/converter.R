
igraph_to_edgelist <- function(g){
    cbind(as_edgelist(g, names = FALSE), E(g)$weight)
}

edgelist_to_igraph <- function(edgelist){
    g <- graph_from_edgelist(edgelist[,1:2], directed=FALSE)
    E(g)$weight <- edgelist[,3]
    return(g)
}



#' Randomizes a weighted graph while keeping the degree distribution constant.
#'
#' Converts the graph to a weighted edge list in NumericMatrix, which is compatible with Rcpp.
#' The Rcpp function "randomize" is called, and then the resulting edge list is converted
#' back into an igraph object.
#' @param g igraph graph
#' @param Q parameter that controls the number of iterations, which will be Q times the order of the graph.
#' @param weight_sel can be either "const_var" or "max_weight".
#' @param lower_bound,upper_bound Bounds to edge weights. The randomization
#' process will avoid steps that would make edge weights fall outside these
#' bounds. Set to NULL for no bound. By default, 0 and NULL respectively.
#' @export
rewireCpp <- function(g, Q=100, weight_sel="const_var", lower_bound=0, upper_bound=NULL){

    
    vertex_names <- igraph::V(g)$name
    edgelist <- igraph_to_edgelist(g)
    if (is.null(upper_bound))
        rewired_edgelist <- randomize(edgelist, Q, weight_sel, lower_bound)
    else
        rewired_edgelist <- randomize(edgelist, Q, weight_sel, lower_bound, upper_bound)
    rewired_g <- edgelist_to_igraph(rewired_edgelist)
    if (!is.null(vertex_names))
        V(rewired_g)$name <- vertex_names
    return ( rewired_g )
}


#' Clustering coefficient of a graph.
#' @param g igraph graph
#' @param upper_bound Upper bound of the edges. 1 by default. If set to NULL, 
#' the maximum edge weight will be taken
#' @export
weighted_clustering_coefficient <- function(g, upper_bound = 1){
    edgelist <- igraph_to_edgelist(g)
    max_weight <- max(E(g)$weight)
    if (is.null(upper_bound) | max_weight > upper_bound) 
        upper_bound <- max_weight
    clustering_coefficient_Rcpp(edgelist, 0, upper_bound)
}

#' Transitivity of a graph.
#' @param g igraph graph
#' @param upper_bound Upper bound of the edges. 1 by default. If set to NULL, 
#' the maximum edge weight will be taken
#' @export
weighted_transitivity <- function(g, upper_bound = 1){
    edgelist <- igraph_to_edgelist(g)
    max_weight <- max(E(g)$weight)
    if (is.null(upper_bound) | max_weight > upper_bound) 
        upper_bound <- max_weight
    transitivity_Rcpp(edgelist, 0, upper_bound)
}


#currently this seems to give the wrong result.
triangle_participation_ratio2 <- function(g){
    edgelist <- igraph_to_edgelist(g)
    triangle_participation_ratio_Rcpp(edgelist)
}




