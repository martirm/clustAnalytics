m_subgraph <- function(G,s){ #calculates the number of edges inside the subgraph induced by the edges in the list s
  ecount(subgraph(G,s))
}

m_subgraph_w <- function(G,s){ #weight of edges inside the subgraph induced by s
  sum(G[s,s])/2
}

degree_w <- function(G,v){ #weighted degree of vertex v
  sum(G[v,])
}

degree_s <- function(G,v,s){ #weighted degree of v over subgraph s (given by a list of vertices)
  sum(G[v,s])
}

median_degree_w <- function(G){
  M <- get.adjacency(G,attr="weight")
  median(apply(M,1,sum))  #computes the median of the degree sequence
}

cs_w <- function(G,s){ #sum of weight of edges connecting vertices of s to the rest of the graph
  sum(G[s,-s]) #sum the elements of rows in s and columns not in s
}


#network metrics
internal_density_s <- function(G,s){
  n <- length(s)
  if (n<=1) return(0)
  else return(2*m_subgraph_w(G,s)/(n*(n-1)))
}

average_degree <- function(G,s){
  2*m_subgraph_w(G,s)/length(s)
}

FOMD_old <- function(G,s,dm="default"){
  if(dm=="default") dm <- median_degree_w(G)  #if median degree has not been given, compute it
  if (length(s)>=2)deg.seq <- apply(G[s,s],1,sum)
  else deg.seq <- 0
  sum(deg.seq>dm)/length(s)
}

expansion <- function(G,s){
  cs_w(G,s)/length(s)
}

cut_ratio <- function(G,s){
  n <- length(s)
  cs_w(G,s)/(n*(gorder(G)-n))
}

conductance <- function(G,s){
  cs_w(G,s)/(2*m_subgraph_w(G,s)+cs_w(G,s))
}

normalized_cut <- function(G,s){
  cs <- cs_w(G,s)
  ms <- m_subgraph_w(G,s)
  m <- sum(M <- get.adjacency(G,attr="weight"))/2
  cs/(2*ms+cs)+cs/(2*(m-ms)+cs)
}

max_odf <- function(G,s){
  first <- TRUE
  M <- 0
  for (i in s){
    out_degree <- sum(G[i,-s]) #out degree (sum of edges leaving s) of vertex i
    odf <- out_degree/degree_w(G,i)
    if (is.na(odf)) next
    if (first | (odf>M)) {
      first <- FALSE
      M <- odf
    }
  }
  return(M)
}

av_odf <- function(G,s){
  sum_odf <- 0
  for (i in s){
    out_degree <- sum(G[i,-s]) #out degree (sum of edges leaving s) of vertex i
    odf <- out_degree/degree_w(G,i)
    sum_odf <- sum_odf + odf
  }
  return(sum_odf/length(s))
}

flake_odf <- function(G,s){
  x <- 0
  for (i in s){ #x will count how many vertices of s have smaller edge weight sum pointing inside than outside
    if ( sum(G[i,s]) < sum(G[i,-s]) ) {
      x <- x+1 
    }
  }
  return(x/length(s))
}

clustering_coef_w <- function(G, n_step=100, w_max=1){
  #delete this in the future. Now replaced by the faster Rcpp implementation.
  mean <-0
  A <- get.adjacency(G,attr="weight")
  if (is.null(w_max)){
      w_max = max(A)
      if (w_max==0)
          return(0)
  }
  for (i in 1:n_step){
    t <- w_max*i/n_step
    At <- A>=t
    Gt <- graph.adjacency(At,mode="undirected")  #unweighted graphs with edges where the weight is above the threshold t
    if (ecount(Gt)>1) {
      trans <- transitivity(Gt)  #if there are no edge we consider the transitivity to be 0 (but igraphs transitivity function would return NaN)
      if (trans!="NaN")mean <- mean+ trans
    }
  }
  return(mean/n_step)
}



clustering_coef_subgraph <- function(G, s, w_max=w_max){
  Gs <- induced_subgraph(G,s)
  weighted_transitivity(Gs, upper_bound=w_max)
}


#' Triangle Participation Ratio (community-wise)
#' 
#' Computes the triangle participation ratio (proportion of vertices that belong 
#' to a triangle). The computation is done to the subgraphs induced by each of 
#' the communities in the given partition.
#' @param g The input graph (as an igraph object). Edge weights and directions are ignored.
#' @param com Community membership vector. Each element corresponds to a vertex
#' of the graph, and contains the index of the community it belongs to.
#' @return A vector containing the triangle participation ratio of each community.
triangle_participation_ratio_communities <- function(g, com){
    TPR_coms_Rcpp(triangles(delete_vertex_attr(g, "name")), com)
}

triangle_participation_ratio <- function(G, s=NULL){
  # computes triangle participation ratio (proportion of vertices that belong to a triangle)
  # computation is done to the subgraph induced by the vertex list s, or to the whole graph if s
  # isn't provided
  # efficiency could be improved by not using the function that counts all triangles
  if (!is.null(s))
    G <- induced_subgraph(G,s)
  a <- adjacent.triangles(G)
  sum(a!=0)/length(a)
}

coverage <- function(G, com){
  # percentage of internal edges respect to total number of edges
  EL <- as_edgelist(G, names=FALSE)
  n_internal <- 0
  f_com <- function(i) com[i]
  EL_com <- apply(EL, c(1,2), f_com)
  internal <- EL_com[,1] == EL_com[,2]
  sum(internal)/length(internal)
}

density_ratio <- function(G, com){
  # new scoring function
  internal_weight <- 0
  for (i in 1:max(com)){
    s <- which(com==i)
    internal_weight <- internal_weight + m_subgraph_w(G, s)
  }
  total_weight <- sum(G[,])/2
  
  sizes <- as.vector(table(com))
  n_total <- sum(sizes)
  potential_internal_n <- function(n) n*(n-1)/2
  potential_external_n <- function(n) n*(n_total-n)
  
  potential_internal <- sum(sapply(sizes, potential_internal_n))
  potential_external <- sum(sapply(sizes, potential_external_n))/2
  
  1 - ((total_weight-internal_weight)/potential_external) / (internal_weight/potential_internal)
}

###################
#other
relabel <- function(c){  
  c2 <- c
  j <- 1
  for (i in unique(c)){
    c2[c==i] <- j
    j <- j+1
  }
  return(c2)
}




##################
weighted_mean_no_NA <- function(v, weights){
    #weighted mean of all the non NA values (the weight corresponding NA values is ignored)
    NA_index <- is.na(v)
    v[NA_index] <- 0
    weights[NA_index] <- 0
    sum(v*weights) /sum(weights)
}

#' Applies function to each subgraph of a graph
#' 
#' @param g igraph graph
#' @param com vector of memberships that determines the subgraphs (i.e. all elements
#' with the same label will form a subgraph).
#' @param f Function to apply. Takes a graph as input and returns a scalar.
#' @return vector with the result of each subgraph
apply_subgraphs <- function(g, com, f, ...){
    labels <- unique(com)  
    f_aux <- function(i) f(induced_subgraph(g, v=which(com==i)), ...)
    return(vapply(labels, f_aux, FUN.VALUE=1))
}


#' #' Cluster scoring functions
#' #' 
#' #' 
#' scoring_functions <- function(G,com,fix_labels=TRUE, globals_only=TRUE, no_clustering_coef=TRUE,
#'                               ignore_NA=TRUE, w_max=1){
#'   if (fix_labels) com <- relabel(com)
#'   n_com <- max(com)
#'   n <- length(V(G))
#'   function_names <- c("internal density","edges inside","av degree","FOMD","expansion",
#'                       "cut ratio","conductance", "norm cut", "max ODF","average ODF","flake ODF","clustering coef","modularity")
#'   D <- data.frame(matrix(nrow=n_com+1,ncol=length(function_names)))
#'   colnames(D) <- function_names
#'   rownames(D) <- c(1:n_com,"global")
#' 
#'   for (i in 1:n_com){
#'     s <- which(com==i) #s contains the indices of vertices belonging to cluster i
#'     if (no_clustering_coef)
#'         cc <- NaN
#'     else
#'         cc <- clustering_coef_subgraph(G, s, w_max=w_max)
#'     D[i,]<- c(internal_density_s(G,s),m_subgraph_w(G,s),average_degree(G,s),FOMD(G,s),expansion(G,s),
#'               cut_ratio(G,s),conductance(G,s),normalized_cut(G,s),max_odf(G,s),av_odf(G,s),flake_odf(G,s),
#'               cc,0)
#'   }
#'   D[n_com+1,]<-colMeans(D[1:n_com,])
#'   weights <- table(com)/n
#'   if (ignore_NA){
#'       D[n_com+1,1:12] <- apply (D[1:n_com,1:12], 2, weighted_mean_no_NA, weights=weights)
#'   }
#'   else{
#'       w_mean <- function(v,weights) sum(v*weights)
#'       D[n_com+1,1:12] <- apply (D[1:n_com,1:12],2,w_mean,weights=weights)
#'   }
#' 
#' 
#'   D["modularity"]<-NaN
#'   D[n_com+1,"modularity"] <- modularity(G,com,weight=E(G)$weight)
#'   #D[n_com+1,"clustering coef"]<- sum(D[1:n_com,"clustering coef"])/ sum(D[1:n_com,"clustering coef"]>0) #mean of non-zero elements
#'   #D[n_com+1,"internal density"]<- sum(D[1:n_com,"internal density"])/ sum(D[1:n_com,"internal density"]>0)
#'   #D[n_com+1,"edges inside"]<- sum(D[1:n_com,"edges inside"])/ sum(D[1:n_com,"edges inside"]>0)
#'   #D[n_com+1,"av degree"]<- sum(D[1:n_com,"av degree"])/ sum(D[1:n_com,"av degree"]>0)
#'   if (globals_only) {
#'       result <- D[n_com+1,]
#'       result["n_communities"] <- length(table(com))
#'       return (result)
#'   }
#'   else return(D)
#' }
#' 
#' scoring_functions_df <- function(G,com, ignore_NA=TRUE, no_clustering_coef=TRUE, type="local", weighted=FALSE){
#'   # type: can be "local" or "global", depending on whether we want a cluster-by cluster or a 
#'   # global analysis
#'   
#'   if (!"weight" %in% names(edge.attributes(G)))
#'     G <- set_edge_attr(G, "weight", value=1)
#'   
#'   n_com <- max(com)
#'   n <- length(V(G))
#'   function_names <- c("size", "internal density","edges inside","av degree","FOMD","expansion",
#'                       "cut ratio","conductance", "norm cut", "max ODF","average ODF","flake ODF", 
#'                       "TPR", "clustering coef","modularity")
#'   D <- data.frame(matrix(nrow=n_com,ncol=length(function_names)))
#'   colnames(D) <- function_names
#'   rownames(D) <- c(1:n_com)
#'   
#'   for (i in 1:n_com){
#'     s <- which(com==i) #s contains the indices of vertices belonging to cluster i
#'     size <- length(s)
#'     Gs <- subgraph(G,s)
#'     if (weighted){
#'       if (no_clustering_coef)
#'         cc <- NaN
#'       else 
#'         cc <- clustering_coef_subgraph(G, s, w_max=w_max)
#'     }
#'     else{
#'       cc <- transitivity(Gs)
#'     }
#'     D[i,]<- c(size, internal_density_s(G,s),m_subgraph_w(G,s),average_degree(G,s),FOMD(G,s),expansion(G,s),
#'               cut_ratio(G,s),conductance(G,s),normalized_cut(G,s),max_odf(G,s),av_odf(G,s),flake_odf(G,s),
#'               triangle_participation_ratio(Gs), cc, NaN)
#'   }
#'   if (type=="local"){
#'     D[,"modularity"] <- modularity(G, com)
#'     D$coverage <- coverage(G, com)
#'     return(D)
#'   }
#'   
#'   
#'   # type=="global" from here
#'   D_glob <- apply(D,2,weighted.mean, w=D$size)
#'   D_glob["size"] <- NaN
#'   D_glob["graph_size"] <- n
#'   D_glob["n_clusters"] <- n_com
#'   D_glob["mean_cluster_size"] <- mean(D$size)
#'   
#'   D_glob["modularity"] <- modularity(G, com)
#'   D_glob["coverage"] <- coverage(G, com)
#'   D_glob["density ratio"] <- density_ratio(G, com)
#'   return(t(as.matrix(D_glob)))
#'   
#'   
#' }

