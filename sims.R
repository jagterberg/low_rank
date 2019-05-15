library(Matrix)
source("helper_functions.R")

n <- 1000
eps <- .3
X <- simulate_noise(n)
l <- spectral_norm(X)
beta = max(abs(X))/l
r <- ceiling(log(n*(n-1))/(eps^2*beta^2))
r <- min(n,r)
Q <- generate_Q(n,r)

k <- 1
while(!is_approximation(X,Q,eps) && k < 50) {
  X <- simulate_noise(n)
  beta = max(abs(X))/spectral_norm(X)
  r <- ceiling(log(n*(n-1))/(eps^2*beta^2))
  r <- min(n,r)
  
  X <- simulate_noise(n)
  Q <- generate_Q(n,r)

  print(paste0("k = ",k))
  k <- k+1
}
