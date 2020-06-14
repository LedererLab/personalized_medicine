```r
######################
## Simulation Study ##
######################

##### Settings

# Load the PAVedr function
source("./PAV.R") # Load PAV.edr function
source("./RidgeCv.R") # Load cross-validation function for Ridge
source("./RidgeEstimator.R") # Load estimation function for Ridge

# Fix a random seed for reproducibility
# seed = 1

# Simulation parameters 
num.runs <- 100 # number of runs of the test (for averaging)
num.obs <- 50 # number of observations
num.par <- 100 # number of parameters
num.testing <- 100 # number of testing data vectors per run

# Generate a set of tuning parameters
num.tuning = 300 # Number of tuning parameters
power <- seq(from = -5, to = 5, length.out = num.tuning)
tuning.parameters <- 10 ^ power

# Choose a test case for the testing data:
# Options: random, real
# random - draw testing vectors by random
# real - Kidney data with 1936 parameters that have p value < 0.05 and random 
# testing vectors
test.case <- "real"

if (test.case == "real"){
  num.runs = 1
}

# Set the signal-to-noise ratio
stn = 0.5

# Set the magnitude of mutual correlations
# If none, then set k = NULL
k = NULL

##### Loop over the number of runs
# Initialize absolute error and computation time arrays
errors.pav  <- matrix(NA, nrow = num.runs, ncol = num.testing)
errors.CV10 <- matrix(NA, nrow = num.runs, ncol = num.testing)
errors.CV5  <- matrix(NA, nrow = num.runs, ncol = num.testing)
errors.CV3  <- matrix(NA, nrow = num.runs, ncol = num.testing)
compTime.pav  <- rep(NA, num.runs)
compTime.CV10 <- rep(NA, num.runs)
compTime.CV5  <- rep(NA, num.runs)
compTime.CV3  <- rep(NA, num.runs)

# Odds Ratio
Odds.ratio <- rep(0, num.runs)


for (run in 1:num.runs) {
  
  if(test.case != "real"){
    cat(sprintf("run %d out of %d runs\n", run, num.runs))
  }
  
  ##### Generate the test setting
  
  cat(paste0("test case is : "), test.case, "\n")
  
  cat(paste0("Magnitude of mutual correlation is : "), k, "\n")
  
  cat("  -> generate test case... ")
  
  if(test.case != "real"){
    
    # design matrix
    if (is.null(k) != 1){
      library(MASS)
      Ones <- rep(1, num.par) %*% t(rep(1, num.par))
      Sigma <- (1 - k) * diag(num.par) + k * Ones
      X <- mvrnorm(n = num.obs, 
                   mu = rep(rnorm(1, 0, 10), num.par), 
                   Sigma = Sigma)
    } else {
      X <- matrix(nrow = num.obs, ncol = num.par)
      for (i in 1:num.obs)
      {
        X[i, ] <- rnorm(num.par, mean = rnorm(1, 0, 10), sd = 1)
      }
    }
    
    # Compute the svd of X
    SVD <- svd(X)
    V <- SVD$v # Right singular vectors
    P <- V %*% t(V) # Projector onto row sapce of X  
    
    # regression vector
    beta <- rnorm(num.par, mean = 0, sd = 1) 
    beta <- as.vector(P %*% beta)
    
    # noise (signal to noise ratio = 0.5)
    u <- rnorm(num.obs, mean = 0, sd = sqrt(stn ^ -1 * var(X %*% beta)))
    
    # outcome
    y <- as.vector(X %*% beta + u)
    
    # Visualization of data
    plot(X %*% beta, y, xlab = "x", ylab = "y")
    
    # Normalize each column of the design matrix
    for (i in 1:num.par)
    {
      nv <- as.vector(X[, i])
      X[, i] <- nv/as.numeric(sqrt(t(nv) %*% nv))
    }
    
  }
  
  ##### Generate test data
  if (test.case == "random") {
    
    ### Random case : 
    Z <- matrix(nrow = num.par, ncol = num.testing)
    for (i in 1:num.par)
    {
      Z[i, ] <- runif(num.testing, -1, 1)
    }
    
  } 
  
  if (test.case == "real") {
    
    # Kidney data example 
    load("./data/E-GEOD-33070_gene.RData")
    load("./data/E-GEOD-33070_clinical.RData")
    load("./data/index_for_parameters.RData")
    
    exp <- as.matrix(exp)
    X <- t(exp[index,])
    num.obs <- dim(X)[1] # Number of observations
    num.par <- dim(X)[2] # Number of regression parameters
    
    # Compute the svd of X
    SVD <- svd(X)
    V <- SVD$v # Right singular vectors
    P <- V %*% t(V) # Projector onto row sapce of X  
    
    # regression vector
    beta <- rnorm(num.par, mean = 0, sd = 1) 
    beta <- as.vector(P %*% beta)
    
    # noise (signal to noise ratio = 0.5)
    u <- rnorm(num.obs, mean = 0, sd = sqrt(stn ^ -1 * var(X %*% beta)))
    
    # outcome
    y <- as.vector(X %*% beta + u)
    
    # Visualization of data
    plot(X %*% beta, y, xlab = "x", ylab = "y")
    
    # Normalize each column of the design matrix
    for (i in 1:num.par)
    {
      nv <- as.vector(X[, i])
      X[, i] <- nv/as.numeric(sqrt(t(nv) %*% nv))
    }
    
    # Generate testing vectors
    Z <- matrix(nrow = num.par, ncol = num.testing)
    for (i in 1:num.par)
    {
      Z[i, ] <- runif(num.testing, -1, 1)
    }
    
  } 
  
  cat("done\n")
  
  ##### Tuning parameter calibration with PAVedr and CV
  
  # Perform the PAVedr pipeline
  cat("  -> compute PAVedr tuning parameters... ")
  # start.time <- Sys.time()
  pav <- PavEdr(y, X, Z, tuning.parameters = tuning.parameters)
  # end.time <- Sys.time()
  compTime.pav[run] <- pav$time  # time in ms
  cat("done\n")
  
  # CV k=10
  cat("  -> compute CV (K=10) tuning parameters... ")
  start.time <- Sys.time()
  CV10.tuning <- RidgeCv(y, X, tuning.parameters = tuning.parameters, K = 10)
  CV10.estimator <- RidgeEstimator(y, X, CV10.tuning)
  end.time <- Sys.time()
  compTime.CV10[run] <- as.numeric(difftime(end.time, start.time, 
  units = "sec"))  # time in ms
  cat("done\n")
  
  # CV k=5
  cat("  -> compute CV (K=5) tuning parameters... ")
  start.time <- Sys.time()
  CV5.tuning <- RidgeCv(y, X, tuning.parameters = tuning.parameters, K = 5)
  CV5.estimator <- RidgeEstimator(y, X, CV5.tuning)
  end.time <- Sys.time()
  compTime.CV5[run] <- as.numeric(difftime(end.time, start.time, 
  units = "sec"))  # time in ms
  cat("done\n")
  
  # CV k=3
  cat("  -> compute CV (K=3) tuning parameters... ")
  start.time <- Sys.time()
  CV3.tuning <- RidgeCv(y, X, tuning.parameters = tuning.parameters, K = 3)
  CV3.estimator <- RidgeEstimator(y, X, CV3.tuning)
  end.time <- Sys.time()
  compTime.CV3[run] <- as.numeric(difftime(end.time, start.time, 
  units = "sec"))  # time in ms
  cat("done\n")
  
  
  ##### Compute the personalized prediction errors 
  cat("  -> compute errors... ")
  for (i in 1:num.testing) {
    z <- as.vector(Z[, i])
    # Absolute error
    errors.pav[run, i] <- abs(t(z) %*% (pav$beta.hat[, i] - beta))
    errors.CV10[run, i] <- abs(t(z) %*% (CV10.estimator - beta)) 
    errors.CV5[run, i] <- abs(t(z) %*% (CV5.estimator - beta))
    errors.CV3[run, i] <- abs(t(z) %*% (CV3.estimator - beta))
  }
  cat("done\n")
  
  if (mean(na.omit(errors.pav)) < mean(na.omit(errors.CV10))){
    Odds.ratio[run] = 1
  }
  
  cat(paste0("PAV \n :", 
         mean(na.omit(errors.pav)),
         "\n",
         "CV10 \n :", 
         mean(na.omit(errors.CV10)), 
         "\n",
         "Odds.ratio \n :",
         mean(Odds.ratio[1:run]),
         "\n"))
  
}  # End loop over the number of runs

##### Data processing

compTime.pav.mean <- mean(compTime.pav)
compTime.pav.scaled <- 1.00
compTime.CV10.scaled <- mean(compTime.CV10) / compTime.pav.mean
compTime.CV5.scaled <- mean(compTime.CV5) / compTime.pav.mean
compTime.CV3.scaled <- mean(compTime.CV3) / compTime.pav.mean


##### Output results

output.Data <- matrix(c(
  mean(errors.pav), mean(errors.CV10), mean(errors.CV5), mean(errors.CV3),
  sd(errors.pav), sd(errors.CV10), sd(errors.CV5), sd(errors.CV3),
  compTime.pav.scaled, compTime.CV10.scaled, compTime.CV5.scaled, 
  compTime.CV3.scaled),
  nrow = 4, ncol = 3)
rownames(output.Data) <- c("Pav.edr", "10-fold CV", "5-fold CV", "3-fold CV")
colnames(output.Data)<- c("mean", "sd", "Scaled-compTime")

cat("Simulation Results : \n")

output.Data
```



