
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
#' @param g \code{igraph} graph, which can be weighted.
#' @param Q Numeric. Parameter that controls the number of iterations, which will be Q times the order of the graph.
#' @param weight_sel can be either "const_var" or "max_weight".
#' @param lower_bound,upper_bound Bounds to the edge weights. The randomization
#' process will avoid steps that would make edge weights fall outside these
#' bounds. Set to NULL for no bound. By default, 0 and NULL respectively.
#' @return The rewired graph.
#' @export
rewireCpp <- function(g, Q=100, weight_sel="max_weight", lower_bound=0, upper_bound=NULL){

    directed <- igraph::is_directed(g)
    vertex_names <- igraph::V(g)$name
    edgelist <- igraph_to_edgelist(g)
    if (is.null(upper_bound))
        rewired_edgelist <- randomize(EdgeList=edgelist, Q=Q, weight_sel=weight_sel, lower_bound=lower_bound, directed=directed)
    else
        rewired_edgelist <- randomize(EdgeList=edgelist, Q=Q, weight_sel=weight_sel, 
                                      lower_bound=lower_bound, upper_bound=upper_bound, directed=directed)
    rewired_g <- edgelist_to_igraph(rewired_edgelist)
    
    # isolated vertices could be lost when creating back the new graph from the edge list
    # we add them back here
    missing_vertices <- gorder(g) - gorder(rewired_g) 
    if (missing_vertices > 0){
        rewired_g <- add_vertices(rewired_g, missing_vertices)
    }
    
    
    if (!is.null(vertex_names))
        V(rewired_g)$name <- vertex_names
    return ( rewired_g )
}


#' Weighted clustering coefficient of a weighted graph.
#' 
#' Weighted clustering Computed using the definition given by McAssey, M. P. and Bijma, F. in 
#' "A clustering coefficient for complete weighted networks" (2015).
#' @param g igraph graph
#' @param upper_bound upper bound to the edge weights used to compute the integral
#' @return The weighted clustering coefficient of the graph (a scalar).
#' @family cluster scoring functions
#' @examples 
#' data(karate, package="igraphdata")
#' weighted_clustering_coefficient(karate)
#' @export
weighted_clustering_coefficient <- function(g, upper_bound=NULL){
    if (gorder(g) < 3){
        return(NA)
    }
    edgelist <- igraph_to_edgelist(g)
    if (is.null(upper_bound))
        upper_bound <- max(E(g)$weight)
    clustering_coefficient_Rcpp(edgelist, 0, upper_bound)
}

#' Weighed transitivity of a weighted graph. 
#' 
#' Computed using the definition given by McAssey, M. P. and Bijma, F. in 
#' "A clustering coefficient for complete weighted networks" (2015).
#' @param g igraph graph
#' @param upper_bound upper bound to the edge weights used to compute the integral
#' @return The weighted transitivity of the graph (a scalar).
#' @family cluster scoring functions
#' @examples 
#' data(karate, package="igraphdata")
#' weighted_transitivity(karate)
#' @export
weighted_transitivity <- function(g, upper_bound=NULL){
    if (gorder(g) < 3){
        return(NA)
    }
    edgelist <- igraph_to_edgelist(g)
    if (is.null(upper_bound))
        upper_bound <- max(E(g)$weight)
    transitivity_Rcpp(edgelist, 0, upper_bound)
}


# Currently this seems to give the wrong result. Until it is fixed, the other 
# implementation will be used
triangle_participation_ratio2 <- function(g){
    edgelist <- igraph_to_edgelist(g)
    triangle_participation_ratio_Rcpp(edgelist)
}




