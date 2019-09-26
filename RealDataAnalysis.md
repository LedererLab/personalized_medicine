```r
###########################
## Kidney data analysis ## 
###########################

##### Settings

# Choose a analysis case for the kidney data:
# Options: in-sample, out-of-sample
# in-sample - In sample prediction for kidney data
# out-of-sample - Out-of-sample leave-one-out prediction for kidney data
analysis.case <- "in-sample"

# Load the PAVedr function
source("./PAV.R") # Load PAV.edr function
source("./RidgeCv.R") # Load cross-validation function for Ridge
source("./RidgeEstimator.R") # Load estimation function for Ridge

# Fix a random seed for reproducibility
# seed = 1

# Load kidney data
load("./data/E-GEOD-33070_gene.RData")
load("./data/E-GEOD-33070_clinical.RData")
exp <- as.matrix(exp)

# Set outcome vector and design matrix
y <- rev(as.vector(weightData$WeightChangePct))
X <- t(exp)

# Normalize design matrix
X.test <- matrix(rep(0,dim(X)[1]*dim(X)[2]),nrow=dim(X)[1],ncol=dim(X)[2])
for (i in 1:dim(X)[2])
{
  nv <- as.vector(X[,i])
  X.test[,i] <- nv/as.numeric(sqrt(t(nv)%*%nv))
}

# Generate a set of tuning parameters
num.tuning = 300 # Number of tuning parameters
power <- seq(from = -5, to = 5, length.out = num.tuning)
tuning.parameters <- 10 ^ power

# Initialize absolute error and computation time arrays
errors.pav  <- rep(NA, dim(X)[1])
errors.CV10 <- rep(NA, dim(X)[1])
errors.CV5  <- rep(NA, dim(X)[1])
errors.CV3  <- rep(NA, dim(X)[1])
compTime.pav  <- rep(NA, dim(X)[1])
compTime.CV10 <- rep(NA, dim(X)[1])
compTime.CV5  <- rep(NA, dim(X)[1])
compTime.CV3  <- rep(NA, dim(X)[1])

if (analysis.case == "in-sample"){
  
  
  
  for (s in 1:dim(X)[1])
  {
    cat(paste0("analysis case is : "), analysis.case, "\n")
    cat(sprintf("run %d out of %d testing vector\n", s, dim(X)[1]))
    
    # Testing vector
    Z <- as.matrix(X.test[s, ])
    test <- as.vector(X.test[s,])
    
    ##### Tuning parameter calibration with PAVedr and CV
    
    # Perform the PAVedr pipeline
    cat("  -> compute PAVedr tuning parameters... ")
    pav <- PavEdr(y, X.test, Z, tuning.parameters = tuning.parameters)
    beta.hat <- pav$beta.hat
    compTime.pav[s] <- pav$time  
    cat("done\n")
    
    # CV k=10
    cat("  -> compute CV (K=10) tuning parameters... ")
    start.time <- Sys.time()
    CV10.tuning <- RidgeCv(y, X.test, tuning.parameters = tuning.parameters, 
    K = 10)
    CV10.estimator <- RidgeEstimator(y, X.test, CV10.tuning)
    end.time <- Sys.time()
    compTime.CV10[s] <- as.numeric(difftime(end.time, start.time, 
    units = "sec"))  
    cat("done\n")
    
    # CV k=5
    cat("  -> compute CV (K=5) tuning parameters... ")
    start.time <- Sys.time()
    CV5.tuning <- RidgeCv(y, X.test, tuning.parameters = tuning.parameters, 
    K = 5)
    CV5.estimator <- RidgeEstimator(y, X.test, CV5.tuning)
    end.time <- Sys.time()
    compTime.CV5[s] <- as.numeric(difftime(end.time, start.time, units = "sec"))  
    cat("done\n")
    
    # CV k=3
    cat("  -> compute CV (K=3) tuning parameters... ")
    start.time <- Sys.time()
    CV3.tuning <- RidgeCv(y, X.test, tuning.parameters = tuning.parameters, 
    K = 3)
    CV3.estimator <- RidgeEstimator(y, X.test, CV3.tuning)
    end.time <- Sys.time()
    compTime.CV3[s] <- as.numeric(difftime(end.time, start.time, units = "sec"))  
    cat("done\n")
    
    errors.pav[s] <- abs(y[s]-t(test)%*%beta.hat)
    errors.CV10[s] <- abs(y[s]-t(test)%*%CV10.estimator)
    errors.CV5[s] <- abs(y[s]-t(test)%*%CV5.estimator)
    errors.CV3[s] <- abs(y[s]-t(test)%*%CV3.estimator)
  
  }
  
} else {
  
  for (s in 1:dim(X)[1])
  {
    cat(paste0("analysis case is : "), analysis.case, "\n")
    cat(sprintf("run %d out of %d testing vector\n", s, dim(X)[1]))
    
    X.run <- X[-s, ]
    y.run <- y[-s]
    
    for (i in 1:dim(X.run)[2])
    {
      nv <- as.vector(X.run[,i])
      X.run[,i] <- nv/as.numeric(sqrt(t(nv)%*%nv))
    }
    
    test <- as.vector(X.test[s,])
    Z <- as.matrix(test)
    
    ##### Tuning parameter calibration with PAVedr and CV
    
    # Perform the PAVedr pipeline
    cat("  -> compute PAVedr tuning parameters... ")
    pav <- PavEdr(y.run, X.run, Z, tuning.parameters = tuning.parameters)
    beta.hat <- as.vector(pav$beta.hat)
    compTime.pav[s] <- pav$time  # time in ms
    cat("done\n")
    
    # CV k=10
    cat("  -> compute CV (K=10) tuning parameters... ")
    start.time <- Sys.time()
    CV10.tuning <- RidgeCv(y.run, X.run, tuning.parameters = tuning.parameters, 
    K = 10)
    CV10.estimator <- RidgeEstimator(y.run, X.run, CV10.tuning)
    end.time <- Sys.time()
    compTime.CV10[s] <- as.numeric(difftime(end.time, start.time, 
    units = "sec"))  
    cat("done\n")
    
    # CV k=5
    cat("  -> compute CV (K=5) tuning parameters... ")
    start.time <- Sys.time()
    CV5.tuning <- RidgeCv(y.run, X.run, tuning.parameters = tuning.parameters, 
    K = 5)
    CV5.estimator <- RidgeEstimator(y.run, X.run, CV5.tuning)
    end.time <- Sys.time()
    compTime.CV5[s] <- as.numeric(difftime(end.time, start.time, units = "sec"))  
    cat("done\n")
    
    # CV k=3
    cat("  -> compute CV (K=3) tuning parameters... ")
    start.time <- Sys.time()
    CV3.tuning <- RidgeCv(y.run, X.run, tuning.parameters = tuning.parameters, 
    K = 3)
    CV3.estimator <- RidgeEstimator(y.run, X.run, CV3.tuning)
    end.time <- Sys.time()
    compTime.CV3[s] <- as.numeric(difftime(end.time, start.time, units = "sec"))  
    cat("done\n")
    
    errors.pav[s] <- abs(y[s]-t(test)%*%beta.hat)
    errors.CV10[s] <- abs(y[s]-t(test)%*%CV10.estimator)
    errors.CV5[s] <- abs(y[s]-t(test)%*%CV5.estimator)
    errors.CV3[s] <- abs(y[s]-t(test)%*%CV3.estimator)
  }
}

##### Data processing

compTime.pav.mean <- mean(as.numeric(compTime.pav))
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



