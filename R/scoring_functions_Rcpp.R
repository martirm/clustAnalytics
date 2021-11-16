# Scoring functions implemented using auxiliary Rcpp code for better performance

#' Auxiliary Functions of a Graph Partition
#'
#' Given a weighted graph and a partition into communities, returns the
#' internal edge weight, the size, and the cut weight for each community.
#' @param g Graph to be analyzed (as an igraph object)
#' @param edelist alternatively, the edgelist of the graph
#' @param com Community membership vector. Each element corresponds to a vertex
#' of the graph, and contains the index of the community it belongs to.
auxiliary_functions <- function(g, com, edgelist){
    if(missing(edgelist)) edgelist <- igraph_to_edgelist(g)
    values <- cluster_auxiliary_values_Rcpp(edgelist, com)
    colnames(values) <- c("m_w", "n", "c_w")
    return(values)
}



#' Maximum, Average, and Flake Out Degree Fractions of a Graph Partition
#'
#' Given a weighted graph and a partition into communities, returns the
#' maximum, average and flake out degree fractions of each community.
#' @param g Graph to be analyzed (as an igraph object)
#' @param edelist alternatively, the edgelist of the graph
#' @param com Community membership vector. Each element corresponds to a vertex
#' of the graph, and contains the index of the community it belongs to.
out_degree_fractions <- function(g, com, edgelist){
    if(missing(edgelist)) edgelist <- igraph_to_edgelist(g)
    odfs <- out_degree_fractions_Rcpp(edgelist, com)
    colnames(odfs) <- c("max ODF","average ODF","flake ODF")
    return(odfs)
}


#' Scoring Functions of a Graph Partition
#'
#' Computes the scoring functions of a graph and its clusters.
#' @param g Graph to be analyzed (as an igraph object). If the edges have a "weight" attribute, those will be used as weights (otherwise, all edges are assumed to be 1).
#' @param com Community membership integer vector. Each element corresponds to a vertex
#' of the graph, and contains the index of the community it belongs to.
#' @param type can be "local" for a cluster by cluster analysis, or "global" for
#' a global analysis of the whole graph partition.
#' @param weighted Is the graph weighted? If it is, doesn't compute TPR score.
#' @param w_max Numeric. Upper bound for edge weights. Should be generally left as default (NULL).
#' Only affects the computation of the clustering coefficient.
#' @param no_clustering_coef  Logical. If TRUE, skips the computation of the clustering
#' coefficient (which can be slow on large graphs).
#'
#' @return If \code{type=="local"}, returns a dataframe with a row for each
#' community, and a column for each score. If \code{type=="global"}, returns a
#' single row with the weighted average scores.
#'
#' @examples
#' data(karate, package="igraphdata")
#' scoring_functions(karate, membership(cluster_louvain(karate)))
#' @export
scoring_functions <- function(g, com, no_clustering_coef=TRUE,
                              type="local", weighted=TRUE, w_max=NULL){

    if (!"weight" %in% names(edge.attributes(g))){
        G <- set_edge_attr(g, "weight", value=1)
        w_max <- 1
    }

    el <- igraph_to_edgelist(g)
    aux_vals <- auxiliary_functions(edgelist = el, com = com)

    n_com <- max(com)
    n <- length(V(g))
    function_names <- c("size", "internal density","edges inside","av degree","FOMD","expansion",
                      "cut ratio","conductance", "norm cut", "max ODF","average ODF","flake ODF",
                      "density ratio", "modularity")
    D <- data.frame(matrix(nrow=n_com,ncol=length(function_names)))
    colnames(D) <- function_names
    rownames(D) <- c(1:n_com)

    # just to make the following lines more readable. (note that _s stands for subgraph)
    n_s <- aux_vals[,"n"]
    m_s <- aux_vals[,"m_w"]
    c_s <- aux_vals[,"c_w"]

    D[,"size"] <- n_s
    D[,"internal density"] <- m_s * 2 / (n_s * (n_s-1))
    D[,"edges inside"] <- m_s
    D[,"av degree"] <- m_s / n_s
    D[,"expansion"] <- c_s / n_s
    D[,"cut ratio"] <- c_s / (n_s * (n - n_s))
    D[,"conductance"] <- c_s / (2*m_s + c_s)
    D[,"norm cut"] <- c_s / (2*m_s + c_s) + c_s / ( 2*(sum(m_s)-m_s) + c_s)
    D[,"density ratio"] <- local_density_ratio_Rcpp(aux_vals)

    D[,c("max ODF","average ODF","flake ODF")] <- out_degree_fractions(edgelist = el, com = com)
    D[,"FOMD"] <- FOMD_Rcpp(el, com)

    if (!weighted){
        D[,"TPR"] <- triangle_participation_ratio_communities(g, com)
    }
    if (!no_clustering_coef) {
        D[,"clustering coef"] <- apply_subgraphs(g, com, weighted_transitivity, upper_bound=w_max)
    }

    internal_weight <- sum(m_s)
    total_weight <- internal_weight + sum(c_s)/2 #total edge weight of the graph (the cut weight is halved because each edge gets counted twice)
    coverage <- internal_weight / total_weight

    if (type=="local"){
        return(D)
    }

    # type=="global" from here
    D_glob <- apply(D,2,weighted.mean, w=D$size, na.rm=TRUE)
    #D_glob["size"] <- NaN
    D_glob["graph_order"] <- n
    D_glob["n_clusters"] <- n_com
    D_glob["mean_cluster_size"] <- mean(D$size)

    D_glob["modularity"] <- modularity(g, com)
    D_glob["coverage"] <- coverage
    D_glob["global density ratio"] <- density_ratio_from_aux(aux_vals)
    return(t(as.matrix(D_glob)))
}


density_ratio_from_aux <- function(aux_vals){
    density_ratio_Rcpp(aux_vals)
}


#' FOMD (Fraction Over Median Degree)
#'
#' Given a weighted graph and a partition into communities, returns the fraction
#' of nodes of each community whose internal degree (i.e. the degree accounting
#' only intra-community edges) is greater than the median degree of the whole
#' graph.
#' @param g Graph to be analyzed (as an igraph object). If the edges have a "weight" 
#' attribute, those will be used as weights.
#' @param edelist alternatively, the edgelist of the graph, as a matrix where the 
#' first two columns to the vertices and the third is the weight of each edge.
#' @param com Community membership integer vector. Each element corresponds to a vertex.
#' @export
FOMD <- function(g, com, edgelist){
    if(missing(edgelist)) edgelist <- igraph_to_edgelist(g)
    FOMD_Rcpp(edgelist, com)
}


#' Internal Density
#' 
#' Internal density of a graph's communities
#' @param g Graph to be analyzed (as an igraph object). If the edges have a "weight" 
#' attribute, those will be used as weights.
#' @param com community membership integer vector. Each element corresponds to a vertex.
#' @return Numeric vector with the internal density of each community.
#' @export
internal_density <- function(g, com, weighted=TRUE, w_max=NULL){
    if (!"weight" %in% names(edge.attributes(g))){
        G <- set_edge_attr(g, "weight", value=1)
        w_max <- 1
    }
    el <- igraph_to_edgelist(g)
    aux_vals <- auxiliary_functions(edgelist = el, com = com)
    # just to make the following lines more readable. (note that _s stands for subgraph)
    n_s <- aux_vals[,"n"]
    m_s <- aux_vals[,"m_w"]
    c_s <- aux_vals[,"c_w"]
    
    return(m_s * 2 / (n_s * (n_s-1)))
}

#' Edges Inside
#' 
#' Number of edges inside a graph's communities, or their accumulated weight if
#' the graph's edges are weighted.
#' @param g Graph to be analyzed (as an igraph object). If the edges have a "weight" 
#' attribute, those will be used as weights.
#' @param com community membership integer vector. Each element corresponds to a vertex.
#' @return Numeric vector with the internal edge weight of each community
#' @export
edges_inside <- function(g, com, weighted=TRUE, w_max=NULL){
    if (!"weight" %in% names(edge.attributes(g))){
        G <- set_edge_attr(g, "weight", value=1)
        w_max <- 1
    }
    el <- igraph_to_edgelist(g)
    aux_vals <- auxiliary_functions(edgelist = el, com = com)

    return(aux_vals[,"m_w"])
}

#' Average Degree
#' 
#' Average degree (weighted degree, if the graph is weighted) of a graph's communities.
#' @param g Graph to be analyzed (as an igraph object). If the edges have a "weight" 
#' attribute, those will be used as weights.
#' @param com community membership integer vector. Each element corresponds to a vertex.
#' @return Numeric vector with the average degree of each community.
#' @export
average_degree <- function(g, com, weighted=TRUE, w_max=NULL){
    if (!"weight" %in% names(edge.attributes(g))){
        G <- set_edge_attr(g, "weight", value=1)
        w_max <- 1
    }
    el <- igraph_to_edgelist(g)
    aux_vals <- auxiliary_functions(edgelist = el, com = com)
    # just to make the following lines more readable. (note that _s stands for subgraph)
    n_s <- aux_vals[,"n"]
    m_s <- aux_vals[,"m_w"]
    return(m_s / n_s)
}

#' Expansion
#' 
#' Given a graph (possibly weighted) split into communities, the expansion of a community
#' is the sum of all edge weights connecting it to the rest of the graph divided by the number
#' of vertices in the community
#' @param g Graph to be analyzed (as an igraph object). If the edges have a "weight" 
#' attribute, those will be used as weights.
#' @param com community membership integer vector. Each element corresponds to a vertex.
#' @return Numeric vector with the expansion of each community.
#' @export
expansion <- function(g, com, weighted=TRUE, w_max=NULL){
    if (!"weight" %in% names(edge.attributes(g))){
        G <- set_edge_attr(g, "weight", value=1)
        w_max <- 1
    }
    el <- igraph_to_edgelist(g)
    aux_vals <- auxiliary_functions(edgelist = el, com = com)
    # just to make the following lines more readable. (note that _s stands for subgraph)
    n_s <- aux_vals[,"n"]
    c_s <- aux_vals[,"c_w"]
    
    return(c_s / n_s)
}

#' Cut Ratio
#' 
#' The cut ratio of a graph's community is the total edge weight connecting the community
#' to the rest of the graph divided by number of unordered pairs of vertices such that one
#' belongs to the community and the other does not.
#' @param g Graph to be analyzed (as an igraph object). If the edges have a "weight" 
#' attribute, those will be used as weights.
#' @param com community membership integer vector. Each element corresponds to a vertex.
#' @return Numeric vector with the cut ratio of each community.
#' @export
cut_ratio <- function(g, com, weighted=TRUE, w_max=NULL){
    if (!"weight" %in% names(edge.attributes(g))){
        G <- set_edge_attr(g, "weight", value=1)
        w_max <- 1
    }
    el <- igraph_to_edgelist(g)
    aux_vals <- auxiliary_functions(edgelist = el, com = com)
    # just to make the following lines more readable. (note that _s stands for subgraph)
    n_s <- aux_vals[,"n"]
    m_s <- aux_vals[,"m_w"]
    c_s <- aux_vals[,"c_w"]
    
    return(c_s / (n_s * (n - n_s)))
}

#' Conductance
#' 
#' Conductance of a graph's communities, which is given by
#' \deqn{\frac{c_s}{2m_s + c_s}},
#' where \eqn{c_s} is the weight of the edges connecting the community s to the rest
#' of the graph, and m_s is the internal weight of the community.
#' 
#' @param g Graph to be analyzed (as an igraph object). If the edges have a "weight" 
#' attribute, those will be used as weights.
#' @param com community membership integer vector. Each element corresponds to a vertex.
#' @return Numeric vector with the conductance of each community.
#' @export
conductance <- function(g, com, weighted=TRUE, w_max=NULL){
    if (!"weight" %in% names(edge.attributes(g))){
        G <- set_edge_attr(g, "weight", value=1)
        w_max <- 1
    }
    el <- igraph_to_edgelist(g)
    aux_vals <- auxiliary_functions(edgelist = el, com = com)
    # just to make the following lines more readable. (note that _s stands for subgraph)
    m_s <- aux_vals[,"m_w"]
    c_s <- aux_vals[,"c_w"]
    
    return(c_s / (2*m_s + c_s))
}

#' Normalized cut
#' 
#' Normalized cut of a graph's communities, which is given by 
#' \deqn{\frac{c_s}{2m_s+c_s}+\frac{c_s}{2(m-m_s)+c_s}},
#' where \eqn{c_s} is the weight of the edges connecting the community s to the rest
#' of the graph, \eqn{m_s} is the internal weight of the community, and \eqn{m} is 
#' the total weight of the network.
#' @param g Graph to be analyzed (as an igraph object). If the edges have a "weight" 
#' attribute, those will be used as weights.
#' @param com community membership integer vector. Each element corresponds to a vertex.
#' @return Numeric vector with the normalized cut of each community.
#' @export
normalized_cut <- function(g, com, weighted=TRUE, w_max=NULL){
    if (!"weight" %in% names(edge.attributes(g))){
        G <- set_edge_attr(g, "weight", value=1)
        w_max <- 1
    }
    el <- igraph_to_edgelist(g)
    aux_vals <- auxiliary_functions(edgelist = el, com = com)
    # just to make the following lines more readable. (note that _s stands for subgraph)
    m_s <- aux_vals[,"m_w"]
    c_s <- aux_vals[,"c_w"]
    
    return( c_s / (2*m_s + c_s) + c_s / ( 2*(sum(m_s)-m_s) + c_s) ) 
}

#' Density Ratio
#' 
#' Internal density of a graph's communities. 
#' @param g Graph to be analyzed (as an igraph object). If the edges have a "weight" 
#' attribute, those will be used as weights.
#' @param com community membership integer vector. Each element corresponds to a vertex.
#' @param type can either be "local" or "global"
#' @return Numeric vector with the internal density of each community.
#' @export
density_ratio <- function(g, com, weighted=TRUE, w_max=NULL, type="local"){
    if (!"weight" %in% names(edge.attributes(g))){
        G <- set_edge_attr(g, "weight", value=1)
        w_max <- 1
    }
    el <- igraph_to_edgelist(g)
    aux_vals <- auxiliary_functions(edgelist = el, com = com)
    if (type=="local"){
        return(local_density_ratio_Rcpp(aux_vals))
    }
    else{
        return(density_ratio_from_aux(aux_vals))
    }
}

#' Max Out Degree Fraction
#'
#' Computes the Maximum Out Degree Fraction (Max ODF) of a graph (which can be weighted) 
#' and its communities.
#' @param g Graph to be analyzed (as an igraph object). If the edges have a "weight" attribute, those will be used as weights (otherwise, all edges are assumed to be 1).
#' @param com Community membership integer vector. Each element corresponds to a vertex
#' of the graph, and contains the index of the community it belongs to.
#' @return Numeric vector with the Max ODF of each community.
#' @export
max_odf <- function(g, com, no_clustering_coef=TRUE,
                              type="local", weighted=TRUE, w_max=NULL){
    
    if (!"weight" %in% names(edge.attributes(g))){
        G <- set_edge_attr(g, "weight", value=1)
        w_max <- 1
    }
    el <- igraph_to_edgelist(g)
    ODFs <- out_degree_fractions(edgelist = el, com = com)
    return(ODFs[,1])
}

#' Average Out Degree Fraction
#'
#' Computes the Average Out Degree Fraction (Average ODF) of a graph (which can be weighted) 
#' and its communities.
#' @param g Graph to be analyzed (as an igraph object). If the edges have a "weight" attribute, those will be used as weights (otherwise, all edges are assumed to be 1).
#' @param com Community membership integer vector. Each element corresponds to a vertex
#' of the graph, and contains the index of the community it belongs to.
#' @return Numeric vector with the Average ODF of each community.
#' @export
average_odf <- function(g, com, no_clustering_coef=TRUE,
                    type="local", weighted=TRUE, w_max=NULL){
    
    if (!"weight" %in% names(edge.attributes(g))){
        G <- set_edge_attr(g, "weight", value=1)
        w_max <- 1
    }
    el <- igraph_to_edgelist(g)
    ODFs <- out_degree_fractions(edgelist = el, com = com)
    return(ODFs[,2])
}

#' Flake Out Degree Fraction
#'
#' Computes the Flake Out Degree Fraction (Max ODF) of a graph (which can be weighted) 
#' and its communities.
#' @param g Graph to be analyzed (as an igraph object). If the edges have a "weight" attribute, those will be used as weights (otherwise, all edges are assumed to be 1).
#' @param com Community membership integer vector. Each element corresponds to a vertex
#' of the graph, and contains the index of the community it belongs to.
#' @return Numeric vector with the Max ODF of each community.
#' @export
max_odf <- function(g, com, no_clustering_coef=TRUE,
                    type="local", weighted=TRUE, w_max=NULL){
    
    if (!"weight" %in% names(edge.attributes(g))){
        G <- set_edge_attr(g, "weight", value=1)
        w_max <- 1
    }
    el <- igraph_to_edgelist(g)
    ODFs <- out_degree_fractions(edgelist = el, com = com)
    return(ODFs[,3])
}
