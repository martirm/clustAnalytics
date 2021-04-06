
#' Test for contingency table approximate counting
#' 
#' Count contingency tables with a singel MC sampling process
MC_test <- function(M){
    H <- 0
    V <- 0
    original_M <- M
    original_M[1,1] <- original_M[1,1] #this is to force original_M to be a deep copy of M. 
    # It would happen automatically when modifying M, except when it's done with Rcpp!!!
    
    v <- rep(0,6)
    
    walk_k_steps(M,0,100)
    while(1){
        dummy <- readline()
        if (dummy == "stop") break
        for(i in 1:100000){
            walk_step(M,0,1)
            if (matrix_equal(M, original_M)) H <- H+1
            V <- V+1
            
            if (M[1,1]==0 & M[1,2]==2) v[1] <- v[1]+1
            else if (M[1,1]==2 & M[2,2]==2) v[2] <- v[2]+1
            else if (M[1,1]==2 & M[2,1]==2) v[3] <- v[3]+1
            else if (M[1,1]==1 & M[2,1]==2) v[4] <- v[4]+1
            else if (M[1,1]==1 & M[3,1]==2) v[5] <- v[5]+1
            else v[6] <- v[6]+1
            
        }
        print(H/V)
        print(v/V)
        
    }
}


matrix_equal <- function(M1, M2){
    for (i in 1:nrow(M1)){
        for (j in 1:ncol(M1)){
            if (M1[i,j] != M2[i,j]) return(FALSE) 
        }
    }
    return(TRUE)
}