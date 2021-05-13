

#' Generates a Barabási-Albert graph with community structure
#' 
#' @param t_max final graph order
#' @param t0 initial
#' @param p vector of label probabilities. If they don't sum 1, they will be scaled accordingly.
#' @param B matrix indicating the affinity of vertices of each label.
#' @param m number of edges added at each step.
#' @param G0 initial graph
#' @param G0_labels labels of the initial graph. If NULL, they will al be set to 1.
#' @param poisson If TRUE, instead of adding m edges incident to each new vertex, 
#' the number of new edges is sampled from a Poisson distribution with parameter m.
#' @param sample_with_replacement If TRUE, allows parallel edges.
#' 
#' @return The resulting graph, as an igraph object. The vertices have a
#' "label" attribute.
#' 
#' 
#' @export
barabasi_albert_blocks <- function(m, p, B, G0=NULL, G0_labels=NULL, t_max,
                                   sample_with_replacement=FALSE,
                                   poisson=FALSE){

    if (is.null(G0)){
        G0 <- generate_G0(n=5*m, p=p, m=m, B=B)
        G0_labels <- V(G0)$label
    }
    G <- G0
    t0 <- gorder(G) + 1
    
    if (poisson){ 
        lambda <- m
    }
    
    new_labels <- sample(1:length(p), t_max-t0+1, replace=TRUE, prob=p)
    if (is.null(G0_labels)){
        G0_labels <- rep(1, gorder(G0))
    }
    labels <- c(G0_labels, new_labels)
    
    degrees <- rep(0, t_max)
    degrees[0:(t0-1)] <- degree(G0)
    
    new_edgelist <- matrix(nrow=(t_max-t0+1)*m , ncol=2)
    i <- 1
    for (t in t0:t_max){
        weights <- degrees[1:(t-1)] * B[labels[1:(t-1)], labels[t]]
        if (poisson){
            m <- rpois(1, lambda)
        }
        new_half_edges <- sample(t-1, size=m, replace=sample_with_replacement, prob=weights)
        
        degrees[t] <- m
        for (j in 1:m){
            v <- new_half_edges[j]
            degrees[v] <- degrees[v] + 1
            new_edgelist[i,] <- c(t, v)
            i <- i+1
        }
    }
    
    G <- add.vertices(G0, nv=t_max-t0+1)
    G <- add.edges(G, as.vector(t(new_edgelist)))
    V(G)$label <- labels
    return(G)
}

generate_G0 <- function(n, p, m, B){
    #TO DO: 
    # -verify that there are enough vertices so that we don't need to add more edges than possible

    
    
    labels <- sample(length(p), n, replace=TRUE, prob=p)
    potential_edges <- matrix(0, nrow=n*(n-1)/2, ncol=3) #first and second columns are adjacent vertices, third is the weight of the edge probability
    k <- 1
    for (i in 1:(n-1)){
        for(j in (i+1):n){
            potential_edges[k,] <- c(i, j, B[labels[i], labels[j]])
            k <- k+1
        }
    }
    selected_edges <- sample(nrow(potential_edges), size=m*n, prob=potential_edges[,3])
    selected_edgelist <- potential_edges[selected_edges,1:2]
    G0 <- graph_from_edgelist(selected_edgelist, directed=FALSE)
    V(G0)$label <- labels
    
    if (any( degree(G0) == 0 )) return(generate_G0(n,p,m,B))
    return(G0)
    
}


