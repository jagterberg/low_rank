#' helper functions for sims.R
#'
#' 
#Rcpp::sourceCpp("matrix_gen.cpp")
library(irlba)
spectral_norm <- function(X) {
  l <- irlba(X,nv =2)
  return(l$d[1])
}

#'  First simulates matrices of the form $A - P$ for  
simulate_noise <- function(n, P = NULL) {
  if (is.null(P)) {
    P <- matrix(rbeta(n^2,2,10000),nrow = n,ncol = n)
    P <- as.matrix(forceSymmetric(P))
  }
  
  A <- generateMatrix(P)

  
  return(A - P)
}


#' Simulate P matrices for the above
#simulate_P_matrix <- function(n) {
  
#}

#' get max norm


#' Compute and check if low-rank approximation works
is_approximation <- function(X,Y,eps,Q = TRUE) {
  if (Q) {
    Y <- get_Y(X,Y)
  }
  
  if(max(abs(X-Y)) <= max(abs(X))*eps) {
    return(TRUE)
  }
  
  return(FALSE)
}


#' generates the n x r matrix of normal random variables
generate_Q <- function(n,r) {
  return( matrix( rnorm(n*r,mean=0,sd=1/(r^(1/2))), n, r) )
}

#' gets the guess for the low-rank approx
get_Y <- function(X,Q) {
  Xsvd <- svd(X)
  new_mat <- Xsvd$u %*% diag(Xsvd$d^(1/2)) %*% Q %*% t(Q) %*% diag(Xsvd$d^(1/2)) %*% Xsvd$v
  return(new_mat)
}


