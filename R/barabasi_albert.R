#' Generates a Barabási-Albert graph with community structure
#'
#' @encoding UTF-8
#' @param t_max maximum value of t (which corresponds to graph order)
#' @param t0 initial t
#' @param p vector of label probabilities. If they don't sum 1, they will be scaled accordingly.
#' @param B matrix indicating the affinity of vertices of each label.
#' @param m number of edges added at each step.
#' @param G0 initial graph
#' @param G0_labels labels of the initial graph. If NULL, they will all be set to 1.
#' @param sample_with_replacement If TRUE, allows parallel edges.
#' @param type Either "Hajek" or "block_first".
#'
#' @return The resulting graph, as an igraph object. The vertices have a
#' "label" attribute.
#'
#'
#' @export
barabasi_albert_blocks <- function(m, p, B, G0=NULL, G0_labels=NULL, t_max,
                                   sample_with_replacement=FALSE,
                                   type="Hajek"){

    if (is.null(G0)){
        G0 <- generate_G0(n=5*m, p=p, m=m, B=B)
        G0_labels <- V(G0)$label
    }
    G <- G0
    t0 <- gorder(G) + 1

    new_labels <- sample(1:length(p), t_max-t0+1, replace=TRUE, prob=p)
    if (is.null(G0_labels)){
        G0_labels <- rep(1, gorder(G0))
    }
    labels <- c(G0_labels, new_labels)

    degrees <- rep(0, t_max)
    degrees[1:(t0-1)] <- degree(G0)

    l <- as.list(t0:t_max)

    if (type == "block_first"){
        fa <-  first_appearance(labels)
    }

    for (t in t0:t_max){
        degrees[t] <- m
        if (type == "Hajek"){
            weights <- degrees[1:(t-1)] * B[labels[1:(t-1)], labels[t]]
            new_half_edges <- sample(t-1, size=m, replace=sample_with_replacement, prob=weights)
            el <- matrix(nrow=m, ncol=2) #edgelist of new edges
            for (j in 1:m){
                v <- new_half_edges[j]
                degrees[v] <- degrees[v] + 1
                el[j,] <- c(t, v)
            }
        }

        if (type == "block_first"){
            populated_communities <- which(fa < t)
            new_half_edge_blocks <- sample(populated_communities, size=m,
                                           replace=TRUE,
                                           prob = B[labels[t], populated_communities])
            el <- matrix(nrow=m, ncol=2) #edgelist of new edges
            k <- 1
            for (j in unique(new_half_edge_blocks)){
                candidates = which(labels[1:t-1]==j)
                n_half_edges = sum(new_half_edge_blocks==j)
                if (length(candidates) == 1){
                    connected_vertices <- rep(candidates, n_half_edges)
                }
                else{
                    connected_vertices <- sample(candidates, size=n_half_edges, prob=degrees[candidates],
                                replace=sample_with_replacement)
                }
                for (v in connected_vertices){
                    degrees[v] <- degrees[v] + 1
                    el[k,] <- c(t, v)
                    k <- k+1
                }
            }
        }
        l[[t-t0+1]] <- el
    }

    complete_edgelist <- do.call(rbind, l)

    G <- add.vertices(G0, nv=t_max-t0+1)
    G <- add.edges(G, as.vector(t(complete_edgelist)))
    V(G)$label <- labels
    return(G)
}

first_appearance <- function(labels){
    # "labels" is a vector of ints with values 1, ..., n
    n <- max(labels)
    first_appearance <- rep(Inf, n)
    for (i in 1:length(labels)){
        l <- labels[i]
        if (is.infinite(first_appearance[l])){
            first_appearance[l] <- i
        }
    }
    return(first_appearance)
}


generate_G0 <- function(n, p, m, B){
    #TO DO:
    # -verify that there are enough vertices so that we don't need to add more edges than possible
    # -add isolated vertices to the graph (if creating from edgelist this doesn't happen)

    n <- round(n)

    labels <- sample(length(p), n, replace=TRUE, prob=p)
    potential_edges <- matrix(0, nrow=n*(n-1)/2, ncol=3) #first and second columns are adjacent vertices, third is the weight of the edge probability
    k <- 1
    for (i in 1:(n-1)){
        for(j in (i+1):n){
            potential_edges[k,] <- c(i, j, B[labels[i], labels[j]])
            k <- k+1
        }
    }
    selected_edges <- sample(nrow(potential_edges), size=round(m*n), prob=potential_edges[,3])
    selected_edgelist <- potential_edges[selected_edges,1:2]
    G0 <- graph_from_edgelist(selected_edgelist, directed=FALSE)
    V(G0)$label <- labels

    if (any( degree(G0) == 0 )) return(generate_G0(n,p,m,B))
    return(G0)

}
