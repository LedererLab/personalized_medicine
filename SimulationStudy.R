######################
## Simulation Study ##
######################

### Settings

# Fix a random seed for reproducibility
# seed = 1

## Simulation parameters 
num.runs <- 3  # number of runs of the test (for averaging)
num.obs <- 50  # number of observations
num.par <- 100  # number of parameters
num.testing <- 100  # number of testing data vectors per run
stn = 0.5  # signal-to-noise ratio
num.tuning = 300 # Number of tuning parameters
# Set the 
# If none, then set k = NULL
k = NULL  # magnitude of mutual correlations

# Choose a test case for the testing data:
# Options: random, real
# random - draw testing vectors by random
# real - Kidney data with 1936 parameters that have p value < 0.05 and random 
# testing vectors
test.case <- "random"

# should the data be plotted?
plot.data <- FALSE


## Load functions and packages
source("./PAV.R") # Load PAV.edr function
source("./RidgeCv.R") # Load cross-validation function for Ridge
source("./RidgeEstimator.R") # Load estimation function for Ridge
source("./Fridge.R") # Load Fridge function
library(MASS)

# Generate a set of tuning parameters
power <- seq(from = -5, to = 5, length.out = num.tuning)
tuning.parameters <- 10 ^ power

# Simulation study based on real data
# [the data for the "random" test case is generated within the run loop.]
if (test.case == "real"){
  # Kidney data example 
  load("./data/E-GEOD-33070_gene.RData")
  load("./data/E-GEOD-33070_clinical.RData")
  load("./data/index_for_parameters.RData")
  
  exp <- as.matrix(exp)
  X <- as.matrix(t(exp[index,]))
  num.obs <- dim(X)[1] # Number of observations
  num.par <- dim(X)[2] # Number of regression parameters
}


##### Loop over the number of runs
# Initialize absolute error and computation time arrays
errors.pav  <- matrix(NA, nrow = num.runs, ncol = num.testing)
errors.CV10 <- matrix(NA, nrow = num.runs, ncol = num.testing)
errors.CV5  <- matrix(NA, nrow = num.runs, ncol = num.testing)
errors.oracle  <- matrix(NA, nrow = num.runs, ncol = num.testing)
errors.fridge  <- matrix(NA, nrow = num.runs, ncol = num.testing)
errors.optimal  <- matrix(NA, nrow = num.runs, ncol = num.testing)

compTime.pav  <- rep(NA, num.runs)
compTime.CV10 <- rep(NA, num.runs)
compTime.CV5  <- rep(NA, num.runs)
compTime.fridge  <- rep(NA, num.runs)

# Odds Ratio
Odds.ratio <- rep(0, num.runs)

# Print simulation setting to screen
cat("Simulation settings: \n",
    "  - test case: ", test.case, "\n",
    "  - n: ", num.obs, ", p: ", num.par, "\n",
    "  - number of testing vectors: ", num.testing, "\n",
    "  - signal-to-noise ratio: ", stn, "\n", sep="")
if(test.case=="random") {
  cat(paste0("  - magnitude of mutual correlation is : "), ifelse(is.null(k), 0, k), "\n")
}
  
# Iterate over number of runs for averaging
for (run in 1:num.runs){
  cat(sprintf("Run %d out of %d runs\n", run, num.runs))
  
  ##### Test case set up
  
  if (test.case == "random"){
    cat("  -> generate test case... ")
    
    # generate design matrix
    if (!is.null(k)){
      # non-zero mutual correlation
      Ones <- rep(1, num.par) %*% t(rep(1, num.par))
      Sigma <- (1 - k) * diag(num.par) + k * Ones
      X <- matrix(nrow = num.obs, ncol = num.par)
      m <- rnorm(num.obs,0,100)
      for (i in 1:num.obs){
        X[i, ] <- mvrnorm(n = 1, 
                          mu = rep(m[i], num.par), 
                          Sigma = Sigma)
      }
    } else {
      # zero mutual correlation
      X <- matrix(nrow = num.obs, ncol = num.par)
      for (i in 1:num.obs)
      {
        X[i, ] <- rnorm(num.par, mean = rnorm(1, 0, 100), sd = 1)
      }
    }
  }
  
  # Normalize each column of the design matrix
  for (i in 1:num.par){
    nv <- as.vector(X[, i])
    X[, i] <- nv/as.numeric(sqrt(t(nv) %*% nv))
  }
  
  # Compute the svd of X
  SVD <- svd(X)
  V <- SVD$v  # Right singular vectors
  P <- V %*% t(V)  # Project onto row space of X  
  
  # regression vector
  beta <- rnorm(num.par, mean = 0, sd = 10) 
  beta <- as.vector(P %*% beta)
  
  # noise (signal to noise ratio = stn)
  u <- rnorm(num.obs, mean = 0, sd = sqrt(stn ^ -1 * var(X %*% beta)))
  
  # outcome
  y <- as.vector(X %*% beta + u)
  
  # Visualization of data
  if(plot.data) plot(X %*% beta, y, xlab = "x", ylab = "y")
  
  
  ## Generate test data
  
  # Random case : 
  Z <- matrix(nrow = num.par, ncol = num.testing)
  for (i in 1:num.testing){
    #Z[i, ] <- runif(num.testing, -1, 1)
    Z[, i] <- rnorm(num.par, 0, 10)
    #Z[, i] <- nv/as.numeric(sqrt(t(nv) %*% nv))
  }
  
  cat("done\n")
  
  
  
  ##### Tuning parameter calibration 
  
  
  ### PAVedr
  
  cat("  -> compute PAVedr tuning parameters... ")
  pav <- PavEdr(y, X, Z, tuning.parameters = tuning.parameters)
  compTime.pav[run] <- pav$time.const  # time in ms
  cat("done\n")
  
  
  ### Oracle and optimal tuning parameters
  
  # Perform the Oracle and Optimal pipeline
  cat("  -> compute Oracle and Optimal tuning parameters... ")
  
  # L2-norm function
  EuclideanNorm <- function(x){
    as.numeric(sqrt(t(x) %*% x))
  }
  
  # Compute the ridge estimator for each tuning parameter
  estimators <- matrix(nrow = num.par, ncol = num.tuning)
  
  for(i in 1:num.tuning) {
    estimators[, i] <- RidgeEstimator(y, X, tuning.parameters[i])
  }
  
  # Compute Optimal tuning parameter
  err <- rep(0, num.testing)
  for (i in 1:num.testing){
    z <- as.vector(Z[, i])
    err[i] <- abs(t(z) %*% (estimators[, i] - beta))
  }
  
  lambda <- rep(0, num.testing)
  for (t in 1:num.testing){
    oracle <- rep(0, length(tuning.parameters))
    bound <- rep(0, length(tuning.parameters))
    bound.optimal <- rep(0, length(tuning.parameters))
    for (o in 1:length(tuning.parameters)){
      z <- as.vector(Z[, t])
      z.norm <- EuclideanNorm(z)
      c <- as.numeric(abs(z %*% estimators[, o]) / 
                        (z.norm * EuclideanNorm(estimators[, o])))
      oracle[o] <- 2 * abs(t(X %*% z) %*% (y - X %*% beta)) / (z.norm * c)
      bound[o] <- c * tuning.parameters[o]
      bound.optimal[o] <- c * oracle[o]
    }
    fit.oracle <- tuning.parameters - oracle
    fit.bound <- bound[fit.oracle >= 0]
    if (isTRUE(length(fit.bound)==0)){
      bound <- bound.optimal
      fit.bound <- bound.optimal
      lambda.index <- min(which(bound == min(fit.bound)))
      lambda[t] <- oracle[lambda.index]
    } else {
      lambda.index <- min(which(bound == min(fit.bound)))
      lambda[t] <- tuning.parameters[lambda.index]
    }
  }
  
  oracle.estimator <- matrix(nrow = num.par, ncol = num.testing)
  for (t in 1:num.testing){
    oracle.estimator[, t] <- RidgeEstimator(y, X, lambda[t])
  }
  cat("done\n")
  
  
  
  ### Fridge
  
  cat("  -> compute Fridge tuning parameters... ")
  start.time <- Sys.time()
  
  Fridge.fit <- fridge(X, y, x0=Z, plug.in='RLOOCV', plot.curve=FALSE)
  fridge.lambda <- Fridge.fit$fridge.tuning
  fridge.estimator <- Fridge.fit$fridge.beta
 
  # fridge.lambda <- rep(0, num.testing)
  # for (t in 1:num.testing){
  #   z <- as.vector(Z[, t])
  #   Fridge.fit <- fridge(X,y,x0=z,plug.in = 'RLOOCV',plot.curve=FALSE)
  #   fridge.lambda[t] <- Fridge.fit$fridge.tuning
  # }
  # 
  # fridge.estimator <- matrix(nrow = num.par, ncol = num.testing)
  # for (t in 1:num.testing){
  #   fridge.estimator[, t] <- RidgeEstimator(y, X, fridge.lambda[t])
  # }
  end.time <- Sys.time()
  compTime.fridge[run] <- Fridge.fit$time.plug.in
  cat("done\n")
  
  
  ### CV (K=10 fold)
  cat("  -> compute CV (K=10) tuning parameters... ")
  start.time <- Sys.time()
  CV10.tuning <- RidgeCv(y, X, tuning.parameters = tuning.parameters, K = 10)
  CV10.estimator <- RidgeEstimator(y, X, CV10.tuning)
  end.time <- Sys.time()
  compTime.CV10[run] <- as.numeric(difftime(end.time, start.time, 
                                            units = "sec"))  # time in ms
  cat("done\n")
  
  ### CV (K=5 fold)
  cat("  -> compute CV (K=5) tuning parameters... ")
  start.time <- Sys.time()
  CV5.tuning <- RidgeCv(y, X, tuning.parameters = tuning.parameters, K = 5)
  CV5.estimator <- RidgeEstimator(y, X, CV5.tuning)
  end.time <- Sys.time()
  compTime.CV5[run] <- as.numeric(difftime(end.time, start.time, 
                                           units = "sec"))  # time in ms
  cat("done\n")
  
  
  
  ##### Personalized prediction errors 
  
  
  cat("  -> compute errors... ")
  for (i in 1:num.testing) {
    z <- as.vector(Z[, i])
    # Absolute error
    errors.pav[run, i] <- abs(t(z) %*% (pav$beta.hat[, i] - beta))
    errors.oracle[run, i] <- abs(t(z) %*% (oracle.estimator[, i] - beta))
    errors.optimal[run, i] <- min(err)
    errors.fridge[run, i] <- abs(t(z) %*% (fridge.estimator[, i] - beta))
    errors.CV10[run, i] <- abs(t(z) %*% (CV10.estimator - beta)) 
    errors.CV5[run, i] <- abs(t(z) %*% (CV5.estimator - beta))
    
  }
  cat("done\n")
  
  if (mean(na.omit(errors.pav)) < mean(na.omit(errors.CV10))){
    Odds.ratio[run] = 1
  }
  
  cat("  -> Mean errors:\n",
      "     - PAV:        ", mean(na.omit(errors.pav)), "\n",
      "     - CV5:        ", mean(na.omit(errors.CV5)), "\n",
      "     - CV10:       ", mean(na.omit(errors.CV10)), "\n",
      "     - Fridge:     ", mean(na.omit(errors.fridge)), "\n",
      "     - Odds.ratio: ",mean(Odds.ratio[1:run]), "\n"
      )
  
}  # End loop over the number of runs




##### Data processing

compTime.pav.mean <- mean(compTime.pav)
compTime.pav.scaled <- 1.00
compTime.fridge.scaled <- mean(compTime.fridge) / compTime.pav.mean
compTime.CV10.scaled <- mean(compTime.CV10) / compTime.pav.mean
compTime.CV5.scaled <- mean(compTime.CV5) / compTime.pav.mean



##### Output results

output.Data <- matrix(c(
  mean(errors.optimal), mean(errors.oracle), mean(errors.pav), mean(errors.fridge), mean(errors.CV5), mean(errors.CV10), 
  sd(errors.optimal), sd(errors.oracle), sd(errors.pav), sd(errors.fridge), sd(errors.CV5), sd(errors.CV10), 
  "None",  "None",  compTime.pav.scaled, compTime.fridge.scaled,  compTime.CV5.scaled, compTime.CV10.scaled),
  nrow = 6, ncol = 3)
rownames(output.Data) <- c("Optimal", "Oracle", "Pav.edr", "Fridge", "5-fold CV", "10-fold CV")
colnames(output.Data)<- c("mean", "sd", "Scaled-compTime")

cat("Simulation Results : \n")
print(output.Data)

output.Data



