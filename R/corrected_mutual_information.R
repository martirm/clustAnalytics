

#' Reduced Mutual Information
#' 
#' Computes the Newman's Reduced Mutual Information (RMI) as defined in 
#' \insertCite{corrected_MI_Newman2020}{clustAnalytics}
#' 
#' The implementation is based on equations 23 (25 for the normalized case) and 29 in (reference here).
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
#' @param alternative_approximation Uses equation 28 instead of 29.
#' 
#' @return The value of Newman's RMI (a scalar).
#' 
#' @section References:
#' \insertRef{corrected_MI_Newman2020}{clustAnalytics}
#' 
#' @export
reduced_mutual_information <- function(c1, c2, base=exp(1), normalized=FALSE,
                                       alternative_approximation=FALSE){
    # alternative_approximation=TRUE uses equation 28 to approximate RMI instead
    # of the approximation in equation 29 (both on Newman's 2020 paper defining 
    # the RMI)
    
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

    if (alternative_approximation){
        v_c_rs <- vector_c_rs(c1, c2)
        sum( lfactorial(v_c_rs) ) / n
        choose2 <- function(k) k * (k-1) / 2
        RMI <- sum( lfactorial(v_c_rs) ) / (n * log(base)) -
               2 / n^3 * sum(sapply(a, choose2)) * sum(sapply(b, choose2)) 
        
    }
    else{
    
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
            warning("error: log_omega<0")
    
        MI <- mutual_information_Cpp(c1, c2, a, b)
        if (base != exp(1)) MI <- MI/log(base)
        RMI <- MI - log_omega/n
    }
    
    if (is.na(MI) | !is.finite(MI)){
        print(MI)
        print("Error with MI")
        print(c1)
        print(c2)
    }
    
    if (is.na(log_omega) | !is.finite(log_omega)){
        print(log_omega)
        print("Error with log_omega")
        print(c1)
        print(c2)
    }
    
    if (!normalized){
        return(RMI)
    }
    else{
        RMI_c1 <- reduced_mutual_information(c1, c1, normalized=FALSE, alternative_approximation=alternative_approximation)
        RMI_c2 <- reduced_mutual_information(c2, c2, normalized=FALSE, alternative_approximation=alternative_approximation)
        
        #cat("MI=", MI, "LO=", log_omega/n, "RMI_c1=",  RMI_c1, "RMI_c2", RMI_c2, "\n")
        #cat("RMI=", 2*(MI-log_omega/n)/(RMI_c1 + RMI_c2), "\n")
        if (2*(RMI)/(RMI_c1 + RMI_c2) > 1){
            save(c1, c2, file="problematic_news_clusterings.RData")
        }
        return(2*(RMI)/(RMI_c1 + RMI_c2))
    }
    
}