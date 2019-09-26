```r
############
## PavEdr ##
############

PavEdr <- function(y, X, Z, tuning.parameters = NULL, num.tuning = 300)
{
  # Description :
  #               Ridge tuning parameter calibration based on personalized 
  #               adaptive validation for Euclidean distance Ridge.
  # Usage : 
  #         PavEdr(y, X, Z, tuning.parameters = NULL, num.tuning = 300)
  # Arguments : 
  #   y : A vector of observations of length n.
  #   X : A design matrix of dimension n * p.
  #   Z : A testing matrix of dimension p * t, where each column is a testing  
  #       vector of dimension p.
  #   tuning.parameters : A vector containing a grid of tuning parameters. The 
  #                       default is set to be NULL. If it is NULL, then a 
  #                       sequence of vector will be generated automatially.
  #   num.tuning : The number of tuning parameters if tuning.parameters = NULL. 
  #                The default value is set to be 300.
  # Returns : 
  #   beta.hat : A matrix of dimension p * t, where each column is the estimator
  #              of regession vector corresponding to the column of Z.
  #   y.hat : A vector of dimension t, where each element is the estimated value 
  #           of Euclidean distance Ridge corresponding to the column of Z.
  #   tuning.chosen : A vector of dimension t, where each element is the optimal 
  #                   tuning parameter chosen by personalized adaptive 
  #                   validation for Euclidean distance Ridge corresponding to 
  #                   the column of Z.
  #   time : The computation time in seconds.
  
  
  
  # Set constants
  num.obs <- as.numeric(dim(X)[1])  # number of observations
  num.par <- as.numeric(dim(X)[2])  # number of parameters
  num.data <- as.numeric(dim(Z)[2])  # number of data vectors
  
  # Utility function
  EuclideanNorm <- function(x) {
    as.numeric(sqrt(t(x) %*% x))
  }
  
  ##### Step 1 : Generate a set of tuning parameters
  
  # Default tuning parameters
  if (is.null(tuning.parameters)) {
    power <- seq(from = -3, to = 3, length.out = num.tuning)
    tuning.parameters <- 10.0^power
  } else {  
    # If tuning.parameters are given set the 
    num.tuning <- length(tuning.parameters)
  }  
  # TODO(YD): Catch errors
  if(num.tuning < 2) {
    stop("number of tuning parameters not sufficient (must be > 1).")
  }
  
  ##### Step 2 : Compute the ridge solution path 
  
  start.time <- Sys.time()
  
  # Singular value decomposition of X
  if (num.obs < num.par){
    SVD <- svd(X)
    U <- SVD$u
    V <- SVD$v
  } else {
    SVD <- svd(X, nu = num.obs, nv = num.par)
    U <- SVD$u
    V <- SVD$v
  }
  
  
  # Compute the ridge estimator for each tuning parameter
  estimators <- matrix(nrow = num.par, ncol = num.tuning)
  
  if (num.obs < num.par){
    D.pseudo.Inverse <- matrix(rep(0, num.obs * num.obs), nrow = num.obs, 
    ncol = num.obs)
  } else {
    D.pseudo.Inverse <- matrix(rep(0, num.obs * num.par), nrow = num.par, 
    ncol = num.obs)
  }
  
  for(i in 1:num.tuning) {
    diag(D.pseudo.Inverse) <- SVD$d / (SVD$d ^ 2 + tuning.parameters[i])
    estimators[, i] <- as.vector(V %*% D.pseudo.Inverse %*% t(U) %*% y)
  }
  
  end.time <- Sys.time()
  
  pav.time.1 <- as.numeric(difftime(end.time, start.time, units = "sec"))  
  
  ##### Step 3 : Transform the ridge tuning parameters to their edr counterparts
  # and sort the tuning parameters
  
  # Transform the ridge tuning parameters to edr tuning parameters
  tuning.edr <- rep(0, num.tuning)
  for (i in 1:num.tuning) {
    tuning.edr[i] <- EuclideanNorm(2 * t(X) %*% (y - X %*% estimators[, i]))
  }
  
  # Initialize output arrays
  tuning.hat <- rep(NA, num.data)
  estimator.hat <- matrix(data = NA, nrow = num.par, ncol = num.data) 
  outcome.hat <- rep(NA, num.data)
  
  
  
  # Loop over all data vectors
  for (k in 1:num.data) {
    
    z <- Z[, k]  # data vector
    z.norm <- EuclideanNorm(z)
    
    as.vector(coef[1:min(num.obs, num.par)]))
    
    start.time <- Sys.time()
    # compute the bounds (1+c[z,r])r/2
    bounds <- rep(0, num.tuning)
    c.bounds <- rep(0, num.tuning)
    for (i in 1:num.tuning) {
      c <- as.numeric(abs(z %*% estimators[, i]) / 
      (z.norm * EuclideanNorm(estimators[, i])))
      bounds[i] <- c * tuning.edr[i]
    }
    
    # compute an ordering such that the bounds increase
    ordering <- sort(bounds, 
                     decreasing   = FALSE, 
                     index.return = TRUE, 
                     method       = "shell")$ix
    bounds.ordered <- bounds[ordering]  # Utility vector of ordered bounds
    estimators.ordered <- estimators[, ordering]
    
    
    
    ##### Step 4 : Use the PAVedr method to compute the optimal tuning parameter
    
    i.optimal <- NaN
    i.optimal.found <- FALSE
    for (i in (num.tuning - 1):1) {
      # Check the binary decision variables for estimators.ordered[, i]
      for (j in (i + 1) : num.tuning) {
        # if it does not hold, then the optimal tuning parameter is 
        # the previous one (i-1)
        if (abs(t(z) %*% (estimators.ordered[, i] - estimators.ordered[, j])) 
            > z.norm * (abs(bounds.ordered[i] + bounds.ordered[j]))) {
          i.optimal <- i+1
          i.optimal.found <- TRUE
          break
        }
      }
      if (i.optimal.found) {
        break
      }
    }
    
    
    
    # If no optimal value has been found according to the PAVedr method
    # choose the smalles tuning parameter
    if(!i.optimal.found) {
      i.optimal <- which.min(tuning.parameters[ordering])
      # TODO(YD): Output warning?
    }
    
    # Output chosen tuning parameter and corresponding estimator and outcome
    tuning.hat[k] <- tuning.parameters[ordering][i.optimal]
    estimator.hat[, k] <- estimators.ordered[, i.optimal]
    outcome.hat[k] <- t(z) %*% estimator.hat[, k]
    
    end.time <- Sys.time()
  }
  
  pav.time.2 <- as.numeric(difftime(end.time, start.time, units = "sec"))  
  
  pav.time <- pav.time.1 + pav.time.2
  
  return(list(beta.hat = estimator.hat, y.hat = outcome.hat, 
              tuning.chosen = tuning.hat, time = pav.time))
}
```
