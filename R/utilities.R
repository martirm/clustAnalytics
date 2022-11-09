
#' Make graph weighted
#' 
#' Given a graph, create a "weight" attribute set to 1 for the edges if it doesn't
#' exist already.
#' @param g igraph graph
#' @return igraph graph with either all edge weights set to 1 (if the original
#' graph was unweighted), or to their original weights if they already existed
#' (in this case, the graph isn't modified at all).
make_graph_weighted <- function(g){
    if (!"weight" %in% names(edge.attributes(g))){
        g <- set_edge_attr(g, "weight", value=1)
        w_max <- 1
    }
    return(g)
}






