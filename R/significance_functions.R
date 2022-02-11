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
#' @param g Graph to be analyzed (as an \code{igraph} object)
#' @param alg_list List of clustering algorithms, which take an \code{igraph} graph as
#' input and return an object of the \code{communities} class.
#' @param w_max Numeric. Upper bound for edge weights. Should be generally left as default (\code{NULL}).
#' @param no_clustering_coef Logical. If \code{TRUE}, skips the computation of the clustering
#' coefficient, which is the most computationally costly of the scoring functions.
#' @param ground_truth Logical. If set to \code{TRUE}, computes the scoring functions for 
#' a ground truth clustering, which has to be provided as \code{gt_clustering}
#' @param gt_clustering Vector of integers that correspond to labels of the ground truth clustering. 
#' Only used if \code{ground_truth} is set to \code{TRUE}.
#' @return
#' A data frame with the values of scoring functions (see \code{\link[clustAnalytics]{scoring_functions}}) 
#' of the clusters obtained by
#' applying the clustering algorithms to the graph.
#' @examples
#' data(karate, package="igraphdata")
#' evaluate_significance(karate)
#' @export
evaluate_significance <- function(g, alg_list=list(Louvain=cluster_louvain, 
                                                   "label prop"= cluster_label_prop, 
                                                   walktrap=cluster_walktrap),
                                  no_clustering_coef=TRUE, ground_truth=FALSE, 
                                  gt_clustering=NULL, w_max=NULL){
    #given an algorithm list and a graph
    c_list <- lapply(alg_list, do.call, list(g)) #clusters graph g by all algorithms in list
    #c_list <- lapply(alg_list, do.call_tryCatch, list(g)) #clusters graph g by all algorithms in list
    memberships_list <- lapply(c_list, membership)
    if (ground_truth){
        memberships_list <- c(memberships_list, list(gt_clustering))
        names(memberships_list)[length(alg_list)+1] <- "ground truth"
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
    return(result)
}



matrix_elementwise_apply <- function (f, A, ...){
    #returns the matrix where element i,j is f(A_ij, B_ij)
    #A and B need to be of the same size.
    M <- matrix(mapply(f, A, ...), ncol=ncol(A))
    rownames(M) <- rownames(A)
    colnames(M) <- colnames(A)
    return(M)
}



merge_percentile <- function(table1, table2, table3){
    # Joins table1 (original scores), table2 (mean rewired scores) and table3
    # percentile of original scores within the corresponding distribution of 
    # rewired scores, with the approppriate suffixes.
    colnames(table2) <- lapply(colnames(table2), paste0, "_r")
    colnames(table3) <- lapply(colnames(table3), paste0, '_percentile')
    cbind(table1, table2, table3)
}

merge4 <- function(table1, table2, table3, digits=3){
    #like merge2, but with percentiles in parentheses
    table1  <- apply(table1, 1:2, format, digits = digits)
    table2  <- apply(table2, 1:2, format, digits = digits)
    table3  <- apply(table3, 1:2, format, digits = digits)
    p <- function(a,b,c) paste0(a, "|", b, "(", c, ")")
    n_row <- nrow(table2)
    n_col <- ncol(table2)
    partial_results <- matrix_elementwise_apply(p, table1[1:n_row, 1:n_col], table2, table3)
    merged_table <- table1
    merged_table[1:n_row, 1:n_col] <- partial_results
    return(merged_table)
}

pull_element <- function(A, i, j) A[i,j]
pull_element_list <- function(l, i,j) sapply(l, pull_element, i=i, j=j)
pull_element_list_nodegen <- function(l, i, j, degen_row="n_clusters"){
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
    ncol_l <- ncol(l[[1]])
    nrow_l <- nrow(l[[1]])
    
    if (ncol(M) > ncol_l){
        M <- M[,1:ncol_l]
    }
    if (nrow(M) >nrow_l){
        M <- M[1:nrow_l,]
    }
    external_conn_scores <- c("expansion", "cut ratio", "conductance", "norm cut", "max ODF", "average ODF", "flake ODF")
    other_scores <- setdiff(rownames(M), external_conn_scores)
    
    element_rank <- function(x, v) dplyr::percent_rank(c(x,v))[1]
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


#' Evaluates the significance of a graph's clusters
#' 
#' Computes community scoring functions to the communities obtained by applying 
#' the given clustering algorithms to a graph. These are compared to the same scores
#' for randomized versions of the graph obtained by a switching algorithm that
#' rewires edges.
#' @inherit evaluate_significance params
#' @param Q Numeric. Parameter that controls the number of iterations of the switching algorithm,
#' which will be Q times the order of the graph.
#' @param lower_bound Numeric. Lower bound to the edge weights. The randomization
#' process will avoid steps that would make edge weights fall outside this
#' bound. It should generally be left as 0 to avoid negative weights.
#' @param weight_sel Can be either \code{const_var} or \code{max_weight}.
#' @param n_reps Number of samples of the rewired graph. 
#' @param ignore_degenerate_cl Logical. If TRUE, when computing the means of the 
#' scoring functions, samples with only one cluster will be ignored.
#' See \link[clustAnalytics]{rewireCpp}.
#' @param w_max Numeric. Upper bound for edge weights. The randomization algorithm will avoid steps that would make 
#' edge weights fall outside this bound. Should be generally left as default (\code{NULL}), unless the network has 
#' by nature or by construction a known upper bound.
#' @param table_style By default returns a table with three columns per algorithm: the original one,
#' the mean of the corresponding rewired scores (suffix "_r") and it's percentile rank within the
#' distribution of rewired scores (suffix "_percentile"). If table_style == "string", instead
#' returns a table with a column per algorithm where each element is of the form "original|rewired(percentile)"
#' @return A matrix with the results of each scoring function and algorithm. See \code{table_style} for details.
#' @export
evaluate_significance_r <- function(g, alg_list=list(Louvain=cluster_louvain, 
                                                     "label prop"= cluster_label_prop, 
                                                     walktrap=cluster_walktrap),
                                    no_clustering_coef=FALSE, 
                                    ground_truth=FALSE, gt_clustering=NULL, 
                                    table_style="default",
                                    ignore_degenerate_cl=TRUE, Q=100, 
                                    lower_bound=0, weight_sel="const_var", 
                                    n_reps=5, w_max=NULL){
    ### weight_sel can be either "const_var" or "max_weight"
    ## ignore_degenerate_cl: if TRUE, when computing the means of the scoring functions, samples with
    ## only one cluster will be ignored.
    upper_bound <- w_max
    table1 <- evaluate_significance(g, alg_list, no_clustering_coef, ground_truth, gt_clustering, 
                                    w_max=w_max)
    if (n_reps == 1){
        rewire_g <- rewireCpp(g, Q, lower_bound=lower_bound, upper_bound=upper_bound, weight_sel=weight_sel)
        table2 <- evaluate_significance(rewire_g, alg_list, no_clustering_coef, ground_truth=FALSE,  
                                        w_max=w_max)
        table_per <- table1
        table_per[,] <- NA
    }
    else{
        # aux <- function(x) rewireCpp(g, Q, lower_bound=lower_bound, 
        #                              upper_bound=upper_bound, weight_sel=weight_sel)
        # rewired_graph_list <- lapply(1:n_reps, aux)
        rewired_graph_list <- replicate(n_reps, rewireCpp(g, Q, lower_bound=lower_bound, 
                                                          upper_bound=upper_bound, 
                                                          weight_sel=weight_sel), 
                                        simplify=FALSE)
        table_list <- lapply(rewired_graph_list, evaluate_significance, alg_list=alg_list,
                             no_clustering_coef=no_clustering_coef, 
                             gt_clustering=FALSE)
        
        mean_table <- Reduce('+', table_list)/n_reps
        # this part is causing errors. Unavailable for now
        # mean_table_nodegen <- matrix_elwise_mean_no_degen(table_list)
        table2 <- mean_table
        # if (ignore_degenerate_cl){
        #     table2[lower_is_better,] <- mean_table_nodegen[lower_is_better,]
        #     table_list <- lapply(table_list, nullify_degen)
        # }
        
        table_per <- percentile_matrix(table1, table_list)
    }
    if (ground_truth){
        table2 <- rbind(table2, NA)
        rownames(table2) <- rownames(table1)
        table_per <- rbind(table_per, NA)
        rownames(table_per) <- rownames(table1)
    }
    
    if (table_style=="string"){ 
        results_table <- merge4(table1, table2, table_per)
        if (ground_truth){
            results_table[,"ground truth"] <- format(table1[,"ground truth"], digits=3)
        }
    }
    else{
        results_table <- merge_percentile(table1, table2, table_per)
    }
    
    return(results_table)
}

higher_is_better <- c("edges inside", "internal density", "av degree", 
                      "FOMD", "clustering coef", "modularity") 
lower_is_better <- c("expansion", "cut ratio", "conductance", "norm cut", 
                     "max ODF", "average ODF", "flake ODF")
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
