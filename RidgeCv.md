```r
###########
## RidgeCv 
###########

# Cross validation for Ridge
RidgeCv <- function(y, X, K, tuning.parameters = NULL, num.tuning = 300)
{
  # Default tuning parameters
  if (is.null(tuning.parameters)) {
    power <- seq(from = -15, to = 150, length.out = num.tuning)
    tuning.parameters <- 10^power
  } else {  
    #If tuning.parameters are given
    num.tuning <- length(tuning.parameters)
  }
  
  n.obs <- as.numeric(dim(X)[1])
  p.var <- as.numeric(dim(X)[2])
  test.number <- floor(n.obs / K)
  
  prediction.error.test <- matrix(nrow = n.obs, ncol = num.tuning)
  
  index <- c(1:n.obs)
  
  #set.seed(seed)
  
  for (k in 1:K)
  {
    if (k == K){
      test.index <- index
    } else {
      test.index <- as.vector(sample(index, test.number))
    }
    index <- index[-which(test.index %in% index)]
    
    train.X <- as.matrix(X[-test.index, ])
    train.y <- as.vector(y[-test.index])
    
    test.X <- as.matrix(X[test.index, ])
    test.y <- as.vector(y[test.index])
    
    if( n.obs < p.var) {
      SVD <- svd(train.X)
      U <- SVD$u
      V <- SVD$v
    } else {
      SVD <- svd(train.X, nu = dim(train.X)[1], nv = p.var)
      U <- SVD$u
      V <- SVD$v
    }
    
    
    # SVD <- svd(train.X)
    # D <- diag(length(SVD$d))
    # diag(D) <- SVD$d
    # U <- SVD$u
    # V <- SVD$v
    
    ridge.estimators <- matrix(nrow = p.var, ncol = num.tuning)
    
    if (dim(train.X)[1] < dim(train.X)[2]){
      D.pseudo.Inverse <- matrix(rep(0, dim(train.X)[1] * dim(train.X)[1]), nrow = dim(train.X)[1], ncol = dim(train.X)[1])
    } else {
      D.pseudo.Inverse <- matrix(rep(0, dim(train.X)[1] * dim(train.X)[2]), nrow = dim(train.X)[2], ncol = dim(train.X)[1])
    }
    
    
    for (tune in 1:num.tuning)
    {
      diag(D.pseudo.Inverse) <- SVD$d / (SVD$d ^ 2 + tuning.parameters[tune])
      ridge.estimators[, tune] <- as.vector(V %*% D.pseudo.Inverse %*% t(U) %*% train.y)
      # D.update <- D / (D ^ 2 + tuning.parameters[tune])
      # ridge.estimators <- cbind(ridge.estimators, 
      #as.vector(V %*% D.update %*% t(U) %*% train.y))
    }
    
    if (k == K) {
      save <- n.obs
    } else {
      save <- (test.number + test.number * (k - 1))
    }
    save.index <- c((1 + test.number * (k - 1)):save)
    prediction.error.test[save.index, ] <- as.matrix(abs(test.y - test.X %*% 
                                                           ridge.estimators))
  }
  prediction.error.sum <- seq(length.out = num.tuning)
  for (s in 1:num.tuning)
  {
    prediction.error.sum[s] <- sum(prediction.error.test[, s])
  }
  
  return(tuning.parameter = tuning.parameters[which.min(prediction.error.sum)])
}
```
