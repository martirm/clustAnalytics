

#' Corrected Mutual Information
#' 
#' Computes the Newman's corrected mutual information (insert Latex reference here)
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
#' 
#' @export
reduced_mutual_information <- function(c1, c2, base=exp(1), normalized=FALSE){
    
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
        #0.5 * log( (gamma(mu*R) * gamma(nu*S)) / ( (gamma(nu)*gamma(R))^S * (gamma(mu)*gamma(S))^R ), base=base)
        0.5 * (lgamma(mu*R) + lgamma(nu*S) - S*(lgamma(nu) + lgamma(R)) - R*(lgamma(mu) + lgamma(S)))

    
    # (this section is replaced by the Rcpp mutual information implementation)
    # c_rs <- matrix(0, nrow=R, ncol=S)
    # 
    # for (i in 1:n){
    #     c_rs[c1[i], c2[i]] <- c_rs[c1[i], c2[i]] + 1
    # }
    # MI <- 0
    # for (i in 1:R){
    #     for (j in 1:S){
    #         P_rs <- c_rs[i,j] / n
    #         if (P_rs != 0){
    #             MI <- MI + P_rs * log(P_rs*n^2 / (a[i]*b[j]), base=base)
    #         }
    #     }
    # }
    
    
    MI <- mutual_information_Cpp(c1, c2, a, b)
    if (base != exp(1)) MI <- MI/log(base)
    
    
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
        return(MI - log_omega/n)
    }
    else{
        return(2*(MI-log_omega/n)/(reduced_mutual_information(c1, c1, normalized=FALSE) + 
                                   reduced_mutual_information(c2, c2, normalized=FALSE)))
    }
    
}