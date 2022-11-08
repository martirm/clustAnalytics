

#' Reduced Mutual Information
#' 
#' Computes the Newman's Reduced Mutual Information (RMI) as defined in 
#' \insertCite{corrected_MI_Newman2020}{clustAnalytics}.
#' 
#' The implementation is based on equations 23 (25 for the normalized case) and 29 in \insertCite{corrected_MI_Newman2020}{clustAnalytics}.
#' The evaluations of the \eqn{\Gamma} functions can get too large and cause overflow 
#' issues in the intermediate steps, so the following term of equation 29:
#' \deqn{\frac{1}{2} \log \frac{\Gamma(\mu R) \Gamma(\nu S)} {(\Gamma(\nu)\Gamma(R))^S (\Gamma(\mu)\Gamma(S))^R } }
#' is rewritten as 
#' \deqn{\frac{1}{2} (\log\Gamma(\mu R) + \log\Gamma(\nu S) - S\log(\Gamma(\nu) - S\log(\Gamma(R) - R\log\Gamma(\mu) - R\log\Gamma(R)  )},
#' and then the function \link[base]{lgamma} is used instead of \link[base]{gamma}.
#' 
#' @param c1,c2 membership vectors
#' @param base base of the logarithms used in the calculations. Changing it only scales the final value. By default set to e=exp(1). 
#' @param normalized If true, computes the normalized version of the corrected mutual information.
#' @param method Can be "hybrid" (default, combines Monte Carlo with analytical formula), "monte_carlo", 
#' approximation1" (appropriate for partitions into many very small clusters), 
#' or "approximation2" (for partitions into few larger clusters).
#' @param warning set to false to ignore the warning.
#' 
#' @return The value of Newman's RMI (a scalar).
#' 
#' @section References:
#' \insertRef{corrected_MI_Newman2020}{clustAnalytics}
#'
#' @export
reduced_mutual_information <- function(c1, c2, base=2, normalized=FALSE,
                                       method="approximation2", warning=TRUE){

    n <- length(c1)
    if (length(c2) != n){
        warning("Membership vectors have different lenghts")
        return()
    }

    a <- count_labels(c1)
    b <- count_labels(c2)
    R <- length(a)
    S <- length(b)

    if (method == "approximation1"){
        v_c_rs <- vector_c_rs(c1, c2)
        sum( lfactorial(v_c_rs) ) / n
        choose2 <- function(k) k * (k-1) / 2
        RMI <- sum( lfactorial(v_c_rs) ) / (n * log(base)) -
               2 / n^3 * sum(sapply(a, choose2)) * sum(sapply(b, choose2)) 
    }
    else{
        MI <- mutual_information_Cpp(c1, c2, a, b)
        if (base != exp(1)) MI <- MI/log(base)
        if (method == "approximation2"){
            log_omega <- log_omega_estimation(c1, c2, base=base)
        }
        else if (method == "hybrid"){
            log_omega <- count_contingency_tables_log(c1, c2, monte_carlo_only=FALSE)
        }
        else if (method == "monte_carlo"){
            log_omega <- count_contingency_tables_log(c1, c2, monte_carlo_only=TRUE)
        }
        else{
            stop("Choose a valid method (either hybrid, approximation1 or approximation2)")
        }
        if (base != exp(1)) log_omega <- log_omega/log(base)
        RMI <- MI - log_omega/n
    }

    
    if (is.na(MI) | !is.finite(MI)){
        warning ("Error with MI. Value: ", MI, call. = TRUE)
    }
    
    if (is.na(log_omega) | !is.finite(log_omega)){
        warning ("Error with log_omega. Value: ", log_omega, call. = TRUE)
    }
    
    if (!normalized){
        return(RMI)
    }
    else{
        RMI_c1 <- reduced_mutual_information(c1, c1, normalized=FALSE, method=method, warning=FALSE, base=base)
        RMI_c2 <- reduced_mutual_information(c2, c2, normalized=FALSE, method=method, warning=FALSE, base=base)
        return(2*(RMI)/(RMI_c1 + RMI_c2))
    }
}

#' Approximation of log(omega(a,b))
#' 
#' Where omega(a,b) is the number of contingency tables with a, b as row and column sums.
#' This approximation is only good for dense tables.
#' @param c1,c2 membership vectors
#' @param base base of the logarithm (e by default)
#' @keywords internal
log_omega_estimation <- function(c1, c2, base=exp(1)){
    n <- length(c1)
    if (length(c2) != n){
        warning("Membership vectors have different lenghts")
        return()
    }
    
    #a <- as.vector(table(c1))
    #b <- as.vector(table(c2))
    a <- count_labels(c1)
    b <- count_labels(c2)
    R <- length(a)
    S <- length(b)
    
    w <- n / (n + 0.5 * R * S)
    x <- (1 - w) / R + (w * a) / n
    y <- (1 - w) / S + (w * b) / n
    mu <- (R + 1) / (R * sum(y^2)) - 1 / R
    nu <- (S + 1) / (S * sum(x^2)) - 1 / S
    
    # Now we can compute the log Omega(a,b) approximation according to equation (29) in the reference
    log_omega <- (R-1) * (S-1) * log(n + 0.5*R*S, base=base) +
        0.5 * (R+nu-2) * sum(log(y, base=base)) +
        0.5 * (S+mu-2) * sum(log(x, base=base)) +
        0.5 * (lgamma(mu*R) + lgamma(nu*S) - S*(lgamma(nu) + lgamma(R)) - R*(lgamma(mu) + lgamma(S))) / log(base)
    if (log_omega <= 0)
        warning("error: log_omega < 0. The analytical approximation is probably not appropriate for this partition.
                An improved implementation of the function will be provided in future versions of the package.")
    
    if (is.na(log_omega) | !is.finite(log_omega)){
        warning ("Error with log_omega. Value: ", log_omega, call. = TRUE)
    }
    return (log_omega)
}

#' Estimates |H_i|/|H_(i+1)| for the first n_rows rows
#' 
#' @param M contingency table
#' @param n_rows number of rows
#' @param error for the convergence of the method
#' @return vector with all the |H_i|/|H_(i+1)| fractions
#' @keywords internal
H_fractions_rows <- function(M, n_rows, error=0.01){
    estimate_H_fractions(M, n_rows, error)
}


#' Natural logarithm of the number of contingency tables
#' 
#' Given a contingency table, returns the natural logarithm of the number of
#' contingency tables that share the same column and row sums. This implementation
#' combines a Markov Chain Monte Carlo approximation with an analytical formula.
#' The input can be either M a contingency table, or two vectors of labels c1 and c2 
#' (in this case, we are counting contingency tables with the same column an row sums as
#' the one produced by c1 and c2)
#' @param c1,c2 membership vectors
#' @param M contingency table
#' @param monte_carlo_only Uses only the Monte Carlo approximation
count_contingency_tables_log <- function(c1, c2, M=NULL, monte_carlo_only=FALSE){
    if (missing(M)){
        M = c_rs_table(c1, c2)
    }
    M <- sort_matrix(M)
    row_sums <- rowSums(M)
    col_sums <- colSums(M)
    p <- c(sum(row_sums==1), sum(col_sums==1))
    if(monte_carlo_only){
        p[1] <- nrow(M)
    }
    
    upper_H_fracs <- H_fractions_rows(M, p[1])
    
    if (p[1]>=nrow(M)-1 | monte_carlo_only){
        ll_H_fracs <- c(1) #ll stands for lower left
        lr_H_fracs_log <- c(0) #lr stands for lower right
    }
    else{
        lM <- M[(p[1]+1):nrow(M),] #lower matrix
        lMt <- t(lM)
        ll_H_fracs <- H_fractions_rows(lMt, p[2])
        
        if(p[2]>=ncol(M)-1){
            lr_H_fracs_log <- c(0)
        }
        else{
            if (missing(c1) & missing(c2)){
                c = contingency_to_membership_vectors(lMt[(p[2]+1):nrow(lMt),] )
                lr_H_fracs_log <- log_omega_estimation(c[[1]], c[[2]])
            }
            else{
                lr_H_fracs_log <- log_omega_estimation(c1, c2)
            }
        }
    }
    return (sum(log(upper_H_fracs), log(ll_H_fracs), lr_H_fracs_log))
}




