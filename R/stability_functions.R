

perturbate_matrix <- function (M, sd=-1, sd_scale=0.1, truncated=FALSE){
    if (sd==-1){
        sd= sqrt(var(as.vector(M))*sd_scale)
    }
    
    n <- dim(M)[1]
    for (i in 1:(n-1)){
        for (j in (i+1):n){
            if (M[i,j] != 0){
                if (truncated) 
                    M[i,j] <- truncnorm::rtruncnorm(1, mean=M[i,j], sd=sd, a=0, b=1)
                else
                    M[i,j] <- truncnorm::rtruncnorm(1, mean=M[i,j], sd=sd, a=0)
                #M[i,j] <- rnorm(1, mean=M[i,j], sd=sd)
                M[j,i] <- M[i,j]
            }
        }
    }
    return (M)
}

cluster_bootstrap_p <- function (returns, n_samples=10, sd=-1, sd_scale=0.1, truncated=FALSE){
    #replace normal distribution with truncated normal (between 0 and 1)
    n <- dim(returns)[2]
    M_dist <- dcor.M(returns[,1:(n/2)])
    gamma.dist=1.48
    c_dist <- HminCpp(M.dist, gamma.dist, itermax=10)
    vi_distances <- vector (mode="numeric", length=n_samples)
    rand_indices <- vi_distances
    for (i in 1:n_samples){
        S_sample <- sample(n/2, n/2, replace=TRUE)
        M_dist_resample <- M_dist[S_sample,S_sample]
        M_p <- perturbate_matrix(M_dist_resample, sd= sd, sd_scale=sd_scale, truncated=truncated)
        c_dist_resample <- HminCpp(M_p, gamma.dist, itermax=10)
        
        c_dist_nd <- c_dist[S_sample] #gets the original clusters for the new sample
        nd <- !duplicated(S_sample) #vector of bools indicating which of the elements are not duplicates
        c_dist_nd <- c_dist_nd[nd]
        c_dist_resample <- c_dist_resample[nd]
        
        vi_distances[i] <- mcclust::vi.dist(c_dist_nd,c_dist_resample)/log(n)
        rand_indices[i] <- fossil::rand.index(c_dist_nd, c_dist_resample)
        
    }
    print(mean(rand_indices))
    return (mean (vi_distances))
}

cluster_bootstrap_igraph <- function (clust_alg, g, n_samples=10, sd=-1, sd_scale=0.1, truncated=FALSE, edgelist=TRUE, perturbate=TRUE, add_self_edges=TRUE){
    n <- gorder(g)
    c_original <- clust_alg(g)
    memb <- membership(c_original)
    
    vi_distances <- vector (mode="numeric", length=n_samples)
    rand_indices <- vi_distances
    for (i in 1:n_samples){
        S_sample <- sample(n, n, replace=TRUE)
        if (edgelist){
            g_resample <- resample_igraph(g, S_sample, perturbate=perturbate, sd=sd, sd_scale=sd_scale, truncated=truncated, add_self_edges=add_self_edges)
        }
        else{
            M <- as_adjacency_matrix(g, type="both", attr = "weight")
            M_resample <- M[S_sample,S_sample]
            M_p <- perturbate_matrix(M_resample, sd= sd, sd_scale=sd_scale, truncated=truncated)
            g_resample <- graph_from_adjacency_matrix(M_p, mode="undirected", weighted=TRUE)
        }
        c_resample <- clust_alg(g_resample)
        
        memb_nd <- memb[S_sample] #gets the original clusters for the new sample
        memb_resample <- membership(c_resample)
        
        vi_distances[i] <- mcclust::vi.dist(memb_nd, memb_resample)/log(n)
        rand_indices[i] <- fossil::rand.index(memb_nd, memb_resample)
        
    }
    results <- c(mean(rand_indices),mean (vi_distances))
    names(results) <- c('Rand Index', 'VI')
    return (c(mean(rand_indices),mean (vi_distances)))
}

jaccard_distances <- function(memb1, memb2, clusters=NULL){
    # computes all jaccard distances between two membership vectors using the criteria described in Hennig2007
    # (for each cluster in memb1, determines its corresponding cluster in memb2 to be the one with smaller 
    # jaccard distance)
    #
    # "clusters" is a vector of indices of the clusters of the output (i.e. which clusters we want to study)
    # this is necessary to keep the outputs consistent when bootstraping (if a cluster doesn't appear in a
    # sample, we want its jaccard value to be 0, not to be ommited)
    
    if (is.null(clusters))
        c1 <- unique(memb1)
    else
        c1 <- clusters
    c2 <- unique(memb2)
    n_union <- function(x1, x2) sum(memb1==x1 | memb2==x2)
    n_intersection <- function(x1, x2) sum(memb1==x1 & memb2==x2)
    unions <- matrix(nrow=length(c1), ncol=length(c2), dimnames=list(c1,c2))
    intersections <- unions
    for (i in 1:length(c1)){
        for (j in 1:length(c2)){
            unions[i,j] <- n_union(c1[i], c2[j])
            intersections[i,j] <- n_intersection(c1[i], c2[j])
        }
    }
    jaccards <- intersections / unions
    return (1-apply(jaccards,1,max))
}

resample_igraph <- function (g, S_sample, perturbate=FALSE, sd=-1, sd_scale=0.1, truncated=FALSE, percentile=0.95, add_self_edges=FALSE){
    n <- gorder(g)
    el <- cbind( get.edgelist(g, names=FALSE) , round( E(g)$weight, 3 ))
    el_s <- resampled_edgelist(el, S_sample)
    
    if(add_self_edges){
        #adds self edges
        a <- sort(el[,3])
        l = length(a)
        a <- a[as.integer(percentile*l):l]
        self_edges <- cbind(1:n,1:n,sample(a,n, replace=TRUE))
        el_s <- rbind(el_s, self_edges)
    }
    
    if (perturbate){
        p <- function(weight){
            if (sd==-1){
                sd= sqrt(var(el[,3])*sd_scale)
            }
            if (truncated) 
                return(truncnorm::rtruncnorm(1, mean=weight, sd=sd, a=0, b=1))
            else
                return(truncnorm::rtruncnorm(1, mean=weight, sd=sd, a=0))
        }
        el_s[,3] <- sapply(el_s[,3], p)
    }
    g_s <- graph_from_edgelist(el_s[,1:2], directed=FALSE)
    E(g_s)$weight <- el_s[,3]
    
    #the bootstrap process could produce some isolated vertices. when calling graph_from_edgelist they are left out of the graph, so we add them back
    missing_vertices = n-gorder(g_s)
    g_s <- add_vertices(g_s, missing_vertices)
    
    return(simplify(g_s))
}

cluster_bootstrap_comparison <- function(g, clust_alg_list, n_samples=10, sd=-1, 
                                         sd_scale=0.1, truncated=FALSE){
    
    results <- sapply(clust_alg_list, cluster_bootstrap_igraph, g=g, 
                      n_samples=n_samples, sd, sd_scale=sd_scale, 
                      truncated=truncated)
    rownames(results) <- c("Rand index", "VI")
    return(results)
}

cluster_statistic <- function(data, sample, memb_original, clust_alg, g, type="global"){
    # clusters the data with clust_alg and returns the VI and rand index of 
    # the results compared to memb_original
    # Intended to use with the boot function
    # type can be either "global" (uses VI and Rand index) or "cluster-wise" (uses Jaccard)
    n <- gorder(g)
    g_resample <- resample_igraph(g, sample)
    g_resample<-set_vertex_attr(g_resample, name="name", value=names(V(g)[sample]))
    memb_nd <- memb_original[sample] #gets the original clusters for the new sample
    #c_resample <- clust_alg(g_resample)
    c_resample <- tryCatch({clust_alg(g_resample)},
                           error=function(e){
                               print("clustering error (algorithm not converging?)")
                               return(NULL)
                           })
    if (is.null(c_resample)) return(c(NaN, NaN, NaN, NaN, NaN))
    memb_resample <- membership(c_resample)
    if (type=="cluster-wise"){
        clusters <- sort(unique(memb_original))
        #clusters <- c(1,2,3) DELETE
        return(jaccard_distances(memb_nd, memb_resample, clusters))
    }
    else{
        statistics <- c(mcclust::vi.dist(memb_nd, memb_resample)/log(n), 
                        clustAnalytics::reduced_mutual_information(memb_nd, memb_resample, 
                                                                   normalized=TRUE), 
                        fossil::rand.index(memb_nd, memb_resample),
                        mclust::adjustedRandIndex(memb_nd, memb_resample),
                        length(unique(memb_resample)))
        names(statistics) <- c("VI", "NRMI", "Rand", "AdRand", "n_clusters")
        return(statistics)
    }
}



boot_cluster <- function(clust_alg, g, R=999, type="global"){
    c_original <- clust_alg(g)
    memb <- membership(c_original)
    vertex_list <- 1:(gorder(g))
    b <- boot::boot(vertex_list, cluster_statistic, R=R, memb_original=memb, 
              clust_alg=clust_alg, g=g, type=type)
}



retry_cluster_leading_eigen <- function(g, itermax=10){
    for (i in 1:itermax){
        try({
            c <- cluster_leading_eigen(g)
            break
            
        })
    }
    return(c)
}

get_statistics <- function(boot_result)    colMeans(boot_result$t, na.rm=TRUE)

statistics_table <- function(b) sapply(b, get_statistics)

#' Performs nonparametric bootstrap to a graph and a list of clustering algorithms
#' 
#' Performs nonparametric bootstrap on a graph's by resampling its vertices and
#' clustering the results using a list of clustering algorithms. 
#' @param alg_list List of igraph clustering algorithms
#' @param g igraph graph object
#' @param R Number of bootstrap replicates.
#' @param type Can be "global" (Variation of Information, Reduced Mutual Information, 
#' and adjusted Rand Index) or "cluster-wise" (Jaccard distance)
#' 
#' @export
boot_alg_list <- function (alg_list = clust_alg_list, g, R=999, return_data=FALSE, type="global"){
    #type can be "global" (VI and Rand) or "cluster-wise" (Jaccard)
    evaluate_boot_g <- function(clust_alg) boot_cluster(clust_alg=clust_alg, g=g, R=R, type=type)
    b <- lapply(alg_list, evaluate_boot_g)
    if(return_data) return(b)
    statistics <- sapply(b,get_statistics)
    rownames(statistics) <- c("VI", "NRMI", "Rand", "AdRand", "n_clusters")
    return(statistics)
}


#' @param t_names names of the statistics in b and b_rnd. Defaults should not be 
#' changed unless the experiments are as well.
#' @param sel_statistics subset of t_names, used to select which statistics we
#' want the table to display.
get_latex_table <- function(b, b_rnd, graph_name = "", caption = NULL, label=NULL,
                            t_names  = c("VI", "RMI", "Rand", "AdRand", "n_clusters"),
                            sel_statistics = c("VI", "RMI", "AdRand", "n_clusters")){
    
    t <- statistics_table(b)
    rownames(t) <- t_names
    t_rnd <- statistics_table(b_rnd)
    rownames(t_rnd) <- sapply(t_names, paste0, "_rnd")
    table <- cbind(t(t), t(t_rnd))
    
    #select desired rows
    all_sel_statistics <- c(sel_statistics, sapply(sel_statistics, paste0, "_rnd"))
    table <- table[,all_sel_statistics]
    
    
    addtorow <- list()
    addtorow$pos <- list(-1,0)
    n_stats <- length(sel_statistics)
    range1 <- paste0("2-", toString(1+n_stats))
    range2 <- paste0(toString(2+n_stats), "-", toString(2*n_stats+1))
    addtorow$command <- c(paste0("& \\multicolumn{",
                                 toString(n_stats),
                                 "}{c}{original} & \\multicolumn{",
                                 toString(n_stats),
                                 "}{c}{randomized} \\\\
                                  \\cmidrule(l){", range1, "} \\cmidrule(l){", range2, "}"),
                          "\\hline\n")
    
    if (is.null(caption)){
        caption <- paste("Mean values of the metrics after bootstrapping with $R=999$ , for both the", 
                         graph_name, "graph and its randomized counterpart, for all tested clustering algorithms")
    }
    if (is.null(label)){
        label <- paste0(graph_name, "_bootstrap")
    }
    
    aux <- paste(replicate(n_stats, "r"), collapse = "")
    align <- paste0("l|", aux, "|", aux)
    l_table <- xtable(table, align = align, caption=caption, label=label)
    print(l_table, add.to.row=addtorow, hline.after=NULL)
}

jaccard_table <- function(b, graph_name = "", caption = NULL, label=NULL){
    if (is.null(caption)){
        caption <- paste("Jaccard distances for each cluster after bootstrapping the", graph_name, "graph with $R=999$")
    }
    if (is.null(label)){
        label <- paste0(graph_name, "_bootstrap_cw")
    }
    
    distances <- function(b_el) colMeans(b_el$t)
    j_list <- lapply(b, distances)
    n_max <- max(unlist(lapply(j_list, length)))
    table <- matrix(nrow=length(b), ncol=n_max, dimnames=list(names(b)))
    for (i in 1:length(j_list)) table[i,1:length(j_list[[i]])] <- j_list[[i]]
    j_table <- xtable(table, caption=caption, label=label)
    print(j_table)
}