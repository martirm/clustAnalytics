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
#' @param g Graph to be analyzed (as an igraph object). It can be a weighted graph, 
#' in which case the weights of the edges
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
#' @param g Graph to be analyzed (as an igraph object)
#' @param edelist alternatively, the edgelist of the graph
#' @param com Community membership vector. Each element corresponds to a vertex
FOMD <- function(g, com, edgelist){
    if(missing(edgelist)) edgelist <- igraph_to_edgelist(g)
    FOMD_Rcpp(edgelist, com)
}



