
#' Computes possible membership vectors from contingency table
#' 
#' Given a contingency table, obtains a possible pair of corresponding labelings.
#' That is, element M[i,j] is the number of elements that belong to community i
#' in the first labeling and j in the second. 
#' @param M the contingency table
#' @return a list containing the two membership vectors
contingency_to_membership_vectors <- function(M){
    l1 <- rep(0, sum(M))
    l2 <- l1
    k <- 1
    for (i in 1:nrow(M)){
        for (j in 1:ncol(M)){
            if (M[i,j]>=1){
                for (l in 1:M[i,j]){
                    l1[k] <- i
                    l2[k] <- j
                    k <- k+1
                }
            }
        }
    }
    return(list(l1,l2))
}


#' Sort matrix
#' 
#' Given a matrix, rearranges rows and columns so that row sums and col sums
#' end up in ascending order.
#' @param M matrix
#' @return rearranged matrix
sort_matrix <- function(M){
    cs <- colSums(M)
    rs <- rowSums(M)
    col_order <- sort(cs, index.return=TRUE)$ix
    row_order <- sort(rs, index.return=TRUE)$ix
    return(M[row_order, col_order])
}



