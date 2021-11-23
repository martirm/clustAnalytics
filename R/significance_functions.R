# library(pbapply)
# library(dplyr)
# library(parallel)
# library(doParallel)




do.call_tryCatch <- function(f, args){
    #same as the do.call, but executes f inside a tryCatch block, and returns NULL 
    #if there is an error 
    g <- args[[1]]
    c <- tryCatch({
        #do.call(f, g)
        return (f(g))
    }, 
    error=function(c){
        print("error")
        message(c)
        return(make_clusters(g, replicate(gorder(g), 1)))
    },
    warning=function(c){
        print("warning")
        message(c)
        return(make_clusters(g, replicate(gorder(g), 1)))
    }
    )
    return(c)
}

#' Evaluates significance of cluster algorithm results on a graph
#' 
#' Given a graph and a list of clustering algorithms, computes several scoring
#' functions on the clusters found by each of the algorithms.
#' @param g Graph to be analyzed (as an igraph object)
#' @param alg_list list of clustering algorithms, which take an \code{igraph} graph as
#' input and return an object of the \code{communities} class.
#' @param no_clustering_coef if \code{TRUE}, skips the computation of the clustering
#' coefficient, which is the most computationally costly of the scoring functions.
#' @param ground_truth if set to\code{TRUE}, computes the scoring functions for 
#' a ground truth clustering, which has to be provided as \code{gt_clustering}
#' @param gt_clustering ground truth clustering. Only used if \code{ground_truth}
#' is set to \code{TRUE}.
#' 
#' @export
evaluate_significance <- function(g, alg_list=clust_alg_list, no_clustering_coef=TRUE, 
                                  ground_truth=FALSE, gt_clustering=NULL, w_max=1){
    #given an algorithm list and a graph
    c_list <- lapply(alg_list, do.call, list(g)) #clusters graph g by all algorithms in list
    #c_list <- lapply(alg_list, do.call_tryCatch, list(g)) #clusters graph g by all algorithms in list
    memberships_list <- lapply(c_list, membership)
    if (ground_truth){
        memberships_list <- c(memberships_list, list(gt_clustering))
        names(memberships_list)[length(clust_alg_list)+1] <- "ground truth"
    }
    # TO DO:
    # change "scoring_functions" to return and use a matrix/vector all along, not a data frame
    # this is currently only a workaround to prevent the conversion of the results to strings by
    # sapply
    aux <- function(mem){
        scores <- scoring_functions(g, mem, no_clustering_coef=no_clustering_coef, w_max=w_max, type="global")
        v <- as.vector(as.matrix(scores))
        names(v) <- colnames(scores)
        if (ground_truth){
            v <- c(v, mcclust::vi.dist(mem, gt_clustering))
            names(v)[length(v)] <- "VIdist_to_GT"
        }
        return(v)
    } 
    result <- sapply(memberships_list, aux)
}



merge1 <- function(table1, table2, digits = 4){
    #option 1 to merge standard and randomized scoring function results.
    #randomized value goes in parentheses under the corresponding original value.
    table1  <- apply(table1, 1:2, format, digits = digits)
    table2  <- apply(table2, 1:2, format, digits = digits)
    parentheses <- function(s) paste0("(",s,")")
    table2 <- apply(table2, 1:2, parentheses)
    full_table <- as.data.frame(matrix(t(cbind(table1,table2)), ncol= ncol(table1), byrow=TRUE))
    full_table <- cbind(c(rbind(rownames(table1), rep("(randomized)",nrow(table1)))), full_table)
    colnames(full_table)[1] <- "score"
    return(full_table)
}


matrix_elementwise_apply <- function (f, A, ...){
    #returns the matrix where element i,j is f(A_ij, B_ij)
    #A and B need to be of the same size.
    M <- matrix(mapply(f, A, ...), ncol=ncol(A))
    rownames(M) <- rownames(A)
    colnames(M) <- colnames(A)
    return(M)
}

merge2 <- function(table1, table2, digits=4){
    table1  <- apply(table1, 1:2, format, digits = digits)
    table2  <- apply(table2, 1:2, format, digits = digits)
    p <- function(a,b) paste0(a, "|", b)
    matrix_elementwise_apply(p, table1, table2)
}

merge3 <- function(table1, table2, suffix="_r"){
    # Joins both tables, adding the suffix to the second matrix colnames 
    # (this could produce duplicates if colnames of table1 already include 
    # the suffix, for now just don't do that)
    # Also includes table1/table2 to get the relative value of each score over
    # its randomized counterpart
    colnames(table2) <- lapply(colnames(table1), paste0, "_r")
    table3 <- table1/table2
    colnames(table3) <- lapply(colnames(table3), paste0, '_rel')
    cbind(table1, table2, table3)
}

merge4 <- function(table1, table2, table3, digits=3){
    #like merge2, but with percentiles in parentheses
    table1  <- apply(table1, 1:2, format, digits = digits)
    table2  <- apply(table2, 1:2, format, digits = digits)
    table3  <- apply(table3, 1:2, format, digits = digits)
    p <- function(a,b,c) paste0(a, "|", b, "(", c, ")")
    matrix_elementwise_apply(p, table1, table2, table3)
}

pull_element <- function(A, i, j) A[i,j]
pull_element_list <- function(l, i,j) sapply(l, pull_element, i=i, j=j)
pull_element_list_nodegen <- function(l, i, j, degen_row="n_communities"){
    # pulls list of elements (i,j) out of each of the matrices of a list, but only those for which
    # ("n_communities", j) is not 1 (degenerate cases)
    l_ij <- pull_element_list(l, i, j)
    l_nj <- pull_element_list(l, degen_row, j)
    l_ij[l_nj != 1]
}




matrix_elwise_mean_no_degen <- function(l){
    #computes the elementwise mean of a list of matrices
    n <- dim(l[[1]])[1]
    m <- dim(l[[1]])[2]
    M <- l[[1]] #we copy the matrix to get the dimensions and dim names, the values will be overwritten
    for (j in 1:m){
        for (i in 1:n){
            l_ij <- pull_element_list_nodegen(l, i, j)
            if (length(l_ij) == 0) 
                M[i,j] <- NaN
            else 
                M[i,j] <- mean(l_ij)
        }
    }
    return(M)
}

nullify_degen <- function(M){
    #turns all elements into NaNs for the columns in which n_cluster is 1
    for (i in 1:ncol(M)){
        if (M["n_communities", i] == 1)
        M[lower_is_better,i] <- NaN
    }
    return(M)
}


percentile_matrix <- function(M, l){
    # from a matrix M and a list of matrices (all of the same dimensions), for
    # each element of M, compute its percentile rank  with respect to all
    # the corresponding values across l.
    
    external_conn_scores <- c("expansion", "cut ratio", "conductance", "norm cut", "max ODF", "average ODF", "flake ODF")
    other_scores <- setdiff(rownames(M), external_conn_scores)
    
    element_rank <- function(x, v) percent_rank(c(x,v))[1]
    M_per <- M
    for (j in 1:(dim(M))[2]){
        for (i in external_conn_scores)
            M_per[i,j] <- element_rank(M[i,j], pull_element_list_nodegen(l, i, j))
    }
    
    for (i in 1:(dim(M)[1])){
        for (j in 1:(dim(M)[2]))
            M_per[i,j]<- element_rank(M[i,j], pull_element_list(l, i, j))
    }
    return(M_per)
    
}


clust_alg_list <- c(cluster_louvain, cluster_label_prop, cluster_walktrap)
names(clust_alg_list) <- c("Louvain", "label prop", "Walktrap")

#' Computes scoring functions of a graph clustering
#' 
#' Computes the scoring functions for both the graph and multiple randomized samples 
#' generated using a switching model
#' @export
evaluate_significance_r <- function(g, alg_list=clust_alg_list, no_clustering_coef=FALSE, 
                                    ground_truth=FALSE, gt_clustering=NULL, table_style=4,
                                    ignore_degenerate_cl=TRUE,
                                    Q=100, lower_bound=0, upper_bound=NULL, weight_sel="const_var", 
                                    n_reps=1, w_max=1, parallel=FALSE){
    ### weight_sel can be either "const_var" or "max_weight"
    ## ignore_degenerate_cl: if TRUE, when computing the means of the scoring functions, samples with
    ## only one cluster will be ignored.
    table1 <- evaluate_significance(g, alg_list, no_clustering_coef, ground_truth, gt_clustering, 
                                    w_max=w_max)
    if (n_reps == 1){
        rewire_g <- rewireCpp(g, Q, lower_bound=lower_bound, upper_bound=upper_bound, weight_sel=weight_sel)
        table2 <- evaluate_significance(rewire_g, alg_list, no_clustering_coef, ground_truth, gt_clustering, 
                                        w_max=w_max)
    }
    else {
        # aux <- function(x) rewireCpp(g, Q, lower_bound=lower_bound, 
        #                              upper_bound=upper_bound, weight_sel=weight_sel)
        # rewired_graph_list <- lapply(1:n_reps, aux)
        rewired_graph_list <- replicate(n_reps, rewireCpp(g, Q, lower_bound=lower_bound, 
                                                          upper_bound=upper_bound, 
                                                          weight_sel=weight_sel), 
                                        simplify=FALSE)
        if (parallel){
            no_cores <- detectCores(logical = TRUE)
            cl <- makeCluster(no_cores-1)  
            registerDoParallel(cl) 
            clusterExport(cl,list('evaluate_significance', 'cluster_leading_eigen'), envir=environment())
            table_list <- parLapply(cl=cl, rewired_graph_list, evaluate_significance, alg_list=alg_list,
                                   no_clustering_coef=no_clustering_coef, ground_truth=ground_truth,
                                   gt_clustering=gt_clustering)
        }
        else{
            table_list <- lapply(rewired_graph_list, evaluate_significance, alg_list=alg_list,
                                   no_clustering_coef=no_clustering_coef, ground_truth=ground_truth,
                                   gt_clustering=gt_clustering)
        }
        
        mean_table <- Reduce('+', table_list)/n_reps
        mean_table_nodegen <- matrix_elwise_mean_no_degen(table_list)
        table2 <- mean_table
        if (ignore_degenerate_cl){
            table2[lower_is_better,] <- mean_table_nodegen[lower_is_better,]
            table_list <- lapply(table_list, nullify_degen)
        }
        
        table_per <- percentile_matrix(table1, table_list)
    }
    
    
    if (table_style==1) results_table <- merge1(table1, table2)
    else if (table_style==2) results_table <- merge2(table1,table2)
    else if (table_style==3) results_table <- merge3(table1, table2)
    else if (table_style==4) results_table <- merge4(table1, table2, table_per)
    else results_table <- table1/table2
    
    if (ground_truth){
        results_table[,"ground truth"] <- format(table1[,"ground truth"], digits=3)
    }
    return(results_table)
}

higher_is_better <- c("edges inside", "internal density", "av degree", "FOMD", "clustering coef", "modularity") 
lower_is_better <- c("expansion", "cut ratio", "conductance", "norm cut", "max ODF", "average ODF", "flake ODF")
add_arrow <- function(s){
    if (s %in% higher_is_better)
        return(paste("\U2191", s))
    else if (s %in% lower_is_better)
        return(paste("\U2193", s))
    else
        return(s)
}
add_arrow_rownames <- function(table){
    rownames(table) <- lapply(rownames(table), add_arrow)
    return(table)
}
