```r
####################
# Ridge estimator ##
####################

RidgeEstimator <- function(y, X, r)
{
  # Description :
  #               The Computation of ridge estimator based on a given tuning
  #               parameter.
  # Usage : 
  #         RidgeEstimato(y, X, r)
  # Arguments : 
  #   y : A vector of observations of length n.
  #   X : A design matrix of dimension n * p.
  #   r : A numerical value representing the tuning parameter for ridge.
  # Returns : 
  #   A vector of dimension p.
  
  ##### Compute the ridge solution by SVD
  
  if (dim(X)[1] < dim(X)[2]){
    SVD <- svd(X)
    U <- SVD$u
    V <- SVD$v
  } else {
    SVD <- svd(X, nu = dim(X)[1], nv = dim(X)[2])
    U <- SVD$u
    V <- SVD$v
  }
  
  if (dim(X)[1] < dim(X)[2]){
    D.pseudo.Inverse <- matrix(rep(0, dim(X)[1] * dim(X)[1]), 
                               nrow = dim(X)[1], 
                               ncol = dim(X)[1])
  } else {
    D.pseudo.Inverse <- matrix(rep(0, dim(X)[1] * dim(X)[2]), 
                               nrow = dim(X)[2], 
                               ncol = dim(X)[1])
  }
  
  
  diag(D.pseudo.Inverse) <- SVD$d / (SVD$d ^ 2 + r)
  return(as.vector(V %*% D.pseudo.Inverse %*% t(U) %*% y))
}
```
