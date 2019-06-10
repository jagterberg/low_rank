library(Rcpp)
library(Matrix)
library(irlba)

cppFunction('NumericMatrix generateMatrix(NumericMatrix P) {
  srand(time(NULL)); 
  int n = P.nrow();
  NumericMatrix A(n,n);
  int br = 0;
            
    for (int i = 0; i < n; i++) {
      for (int j = 0; j <= i ; j++) {
        double r = ((double) rand() / (RAND_MAX));
          if (r >= 1 - P(i,j)) {
              br = 1;
          }
          A(i,j) = br;
          A(j,i) = A(i,j);
          br = 0;
      }
    }
    
  return A;
}
')

#' Can give either X or Xsvd
spectral_norm <- function(X = NULL,Xsvd = NULL) { 
  if (is.null(Xsvd)) {
    l <- irlba(X,nv =2)
    return(l$d[1])
  } else {
    return(Xsvd$d[1])
  }
  
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

#' Compute and check if low-rank approximation works
#' Give either X or Xsvd directly
#' Give either Y or give Q so that Y is computed from Q
is_approximation <- function(X = NULL,Y = NULL,eps,Q = NULL,Xsvd = NULL) {
  # most efficient:
  if (!is.null(Q) & !is.null(Xsvd)) {
    Y <- get_Y(Xsvd=Xsvd,Q=Q)
  #less efficient (needs svd of X)
  } else if (is.null(Q) & !is.null(Xsvd)) {
    Y <- get_Y(Y=Y,Xsvd = Xsvd)
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
get_Y <- function(X = NULL,Q,Xsvd=NULL) {
  if (is.null(Xsvd)) {
    Xsvd <- svd(X)
  }
  new_mat <- Xsvd$u %*% diag(Xsvd$d^(1/2)) %*% Q %*% t(Q) %*% diag(Xsvd$d^(1/2)) %*% Xsvd$v
  return(new_mat)
}


n <- 10000
eps <- .1
kmax <- 200


print("Simulating a good X matrix")
X <- simulate_noise(n)
Xsvd <- svd(X)
l <- spectral_norm(Xsvd=Xsvd)
beta = max(abs(X))/l
r <- ceiling(log(n*(n-1))/(eps^2*beta^2))
r <- min(n,r)
j <- 1
jmax <- 200
while(r == n && j < jmax) {
  X <- simulate_noise(n)
  Xsvd <- svd(X)
  l <- spectral_norm(Xsvd=Xsvd)
  beta = max(abs(X))/l
  r <- ceiling(log(n*(n-1))/(eps^2*beta^2))
  r <- min(n,r)
  j <- j+1
  print(paste0("Trying new X for the ",j, "th time of ",jmax))
  
}

print("Trying a Q matrix now...")
Q <- generate_Q(n,r)
k <- 1
while(!is_approximation(X=X,Xsvd=Xsvd,Q=Q,eps=eps) && k < kmax) {
  #X <- simulate_noise(n)
  #Xsvd <- svd(X)
  #l <- spectral_norm(Xsvd=Xsvd)
  #beta = max(abs(X))/l
  #r <- ceiling(log(n*(n-1))/(eps^2*beta^2))
  #r <- min(n,r)
  Q <- generate_Q(n,r)
  
  print(paste0("k = ",k, " out of ", kmax))
  k <- k+1
}

if (k == kmax) {
  final_val <- NULL
  print("simulation unsuccessful")
} else {
  final_val <- list(X,Q,r,eps)
  save(final_val,file=paste0("n",n,",r",r,",eps",eps,".RData"))
  print("saved final_val")
}









