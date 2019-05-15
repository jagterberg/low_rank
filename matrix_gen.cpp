#include <Rcpp.h>
using namespace Rcpp;

#include <time.h>
#include <stdlib.h>

  // Initialization, should only be called once.

// [[Rcpp::export]]
NumericMatrix generateMatrix(NumericMatrix P) {
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
