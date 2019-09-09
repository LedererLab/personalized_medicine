##################
# Ridge estimator
##################

RidgeEstimator <- function(y, X, r)
{
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
