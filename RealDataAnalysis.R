###########################
## Kidney data analysis ## 
###########################

##### Settings

# Choose a analysis case for the kidney data:
# Options: in-sample, out-of-sample
# in-sample - In sample prediction for kidney data
# out-of-sample - Out-of-sample leave-one-out prediction for kidney data
analysis.case <- "out-of-sample"

# Load the PAVedr function
source("./PAV.R") # Load PAV.edr function
source("./RidgeCv.R") # Load cross-validation function for Ridge
source("./RidgeEstimator.R") # Load estimation function for Ridge
source("./Fridge.R") # Load Fridge function

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
errors.fridge  <- rep(NA, dim(X)[1])

compTime.pav  <- rep(NA, dim(X)[1])
compTime.CV10 <- rep(NA, dim(X)[1])
compTime.CV5  <- rep(NA, dim(X)[1])
compTime.fridge <- rep(NA, dim(X)[1])

cat(paste0("analysis case is : "), analysis.case, "\n")
if (analysis.case == "in-sample"){
  ##### in-sample prediction
  
  for (s in 1:dim(X)[1]){
    cat(sprintf("run %d out of %d testing vector\n", s, dim(X)[1]))
    
    # Testing vector
    Z <- as.matrix(X.test[s, ])
    test <- as.vector(X.test[s,])
    
    
    ##### Tuning parameter calibration
    
    ### PAVedr
    cat("  -> compute PAVedr tuning parameters... ")
    start.time <- Sys.time()
    pav <- PavEdr(y, X.test, Z, tuning.parameters = tuning.parameters)
    beta.hat <- pav$beta.hat
    compTime.pav[s] <- as.numeric(difftime(Sys.time(), start.time, 
                                           units = "sec"))
    cat("done\n")
    
    ### Fridge
    cat("  -> compute Fridge tuning parameters... ")
    start.time <- Sys.time()
  
    Fridge.fit <- fridge(X,y,x0=Z,plug.in = 'RLOOCV',plot.curve=FALSE)
    fridge.lambda <- Fridge.fit$fridge.tuning
    
    fridge.estimator <- RidgeEstimator(y, X.test, fridge.lambda)
    
    compTime.fridge[s] <- as.numeric(difftime(Sys.time(), start.time, 
                                              units = "sec"))
    cat("done\n")
    
    
    ### CV (K=10 fold)
    cat("  -> compute CV (K=10) tuning parameters... ")
    start.time <- Sys.time()
    CV10.tuning <- RidgeCv(y, X.test, tuning.parameters = tuning.parameters, 
                           K = 10)
    CV10.estimator <- RidgeEstimator(y, X.test, CV10.tuning)
    compTime.CV10[s] <- as.numeric(difftime(Sys.time(), start.time, 
                                            units = "sec"))  
    cat("done\n")
    
    ### CV (K=5 fold)
    cat("  -> compute CV (K=5) tuning parameters... ")
    start.time <- Sys.time()
    CV5.tuning <- RidgeCv(y, X.test, tuning.parameters = tuning.parameters, 
                          K = 5)
    CV5.estimator <- RidgeEstimator(y, X.test, CV5.tuning)
    compTime.CV5[s] <- as.numeric(difftime(Sys.time(), start.time, 
                                           units = "sec"))  
    cat("done\n")
    
    errors.fridge[s] <- abs(y[s]-t(test)%*%fridge.estimator)
    errors.pav[s] <- abs(y[s]-t(test)%*%beta.hat)
    errors.CV10[s] <- abs(y[s]-t(test)%*%CV10.estimator)
    errors.CV5[s] <- abs(y[s]-t(test)%*%CV5.estimator)
    
    cat("  -> Mean errors:\n",
        "     - PAV:        ", mean(na.omit(errors.pav)), "\n",
        "     - CV5:        ", mean(na.omit(errors.CV5)), "\n",
        "     - CV10:       ", mean(na.omit(errors.CV10)), "\n",
        "     - Fridge:     ", mean(na.omit(errors.fridge)), "\n"
    )
  }
  # end of in-sample prediction
} else {
  
  ##### Out-of-sample prediction
  
  for (s in 1:dim(X)[1])
  {
    cat(sprintf("run %d out of %d testing vector\n", s, dim(X)[1]))
    
    X.run <- X.test[-s, ]
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
    start.time <- Sys.time()
    pav <- PavEdr(y.run, X.run, Z, tuning.parameters = tuning.parameters)
    beta.hat <- as.vector(pav$beta.hat)
    compTime.pav[s] <- as.numeric(difftime(Sys.time(), start.time, 
                                           units = "sec"))  
    cat("done\n")
    
    # Fridge
    cat("  -> compute Fridge tuning parameters... ")
    start.time <- Sys.time()
    
    Fridge.fit <- fridge(X,y,x0=Z,plug.in = 'RLOOCV',plot.curve=FALSE)
    fridge.lambda <- Fridge.fit$fridge.tuning
    
    fridge.estimator <- RidgeEstimator(y.run, X.run, fridge.lambda)
    
    compTime.fridge[s] <- as.numeric(difftime(Sys.time(), start.time, 
                                              units = "sec"))
    cat("done\n")
    
    # CV k=10
    cat("  -> compute CV (K=10) tuning parameters... ")
    start.time <- Sys.time()
    CV10.tuning <- RidgeCv(y.run, X.run, tuning.parameters = tuning.parameters, 
                           K = 10)
    CV10.estimator <- RidgeEstimator(y.run, X.run, CV10.tuning)
    compTime.CV10[s] <- as.numeric(difftime(Sys.time(), start.time, 
                                            units = "sec"))  
    cat("done\n")
    
    # CV k=5
    cat("  -> compute CV (K=5) tuning parameters... ")
    start.time <- Sys.time()
    CV5.tuning <- RidgeCv(y.run, X.run, tuning.parameters = tuning.parameters, 
                          K = 5)
    CV5.estimator <- RidgeEstimator(y.run, X.run, CV5.tuning)
    compTime.CV5[s] <- as.numeric(difftime(Sys.time(), start.time, units = "sec"))  
    cat("done\n")
    
    errors.fridge[s] <- abs(y[s]-t(test)%*%fridge.estimator)
    errors.pav[s] <- abs(y[s]-t(test)%*%beta.hat)
    errors.CV10[s] <- abs(y[s]-t(test)%*%CV10.estimator)
    errors.CV5[s] <- abs(y[s]-t(test)%*%CV5.estimator)
    
    cat("  -> Mean errors:\n",
        "     - PAV:        ", mean(na.omit(errors.pav)), "\n",
        "     - CV5:        ", mean(na.omit(errors.CV5)), "\n",
        "     - CV10:       ", mean(na.omit(errors.CV10)), "\n",
        "     - Fridge:     ", mean(na.omit(errors.fridge)), "\n"
        )
    
  }
}

##### Data processing

compTime.pav.mean <- mean(as.numeric(compTime.pav))
compTime.pav.scaled <- 1.00
compTime.CV10.scaled <- mean(compTime.CV10) / compTime.pav.mean
compTime.CV5.scaled <- mean(compTime.CV5) / compTime.pav.mean
compTime.fridge.scaled <- mean(compTime.fridge) / compTime.pav.mean


##### Output results

output.Data <- matrix(c(
  mean(errors.pav), mean(errors.fridge), mean(errors.CV10), mean(errors.CV5),
  sd(errors.pav), sd(errors.fridge), sd(errors.CV10), sd(errors.CV5),
  compTime.pav.scaled, compTime.fridge.scaled, compTime.CV10.scaled, compTime.CV5.scaled),
  nrow = 4, ncol = 3)
rownames(output.Data) <- c("Pav.edr", "Fridge", "10-fold CV", "5-fold CV")
colnames(output.Data)<- c("mean", "sd", "Scaled-compTime")

results <- data.frame(
  "mean"=c(mean(errors.pav), mean(errors.fridge), mean(errors.CV5), mean(errors.CV10)),
  "sd"=c(sd(errors.pav), sd(errors.fridge), sd(errors.CV5), sd(errors.CV10)), 
  "comp.time"=c(mean(compTime.pav), mean(compTime.fridge), mean(compTime.CV10), mean(compTime.CV5)),
  row.names=c("Pav.edr", "Fridge", "5-fold CV", "10-fold CV")
)

cat("Simulation Results : \n")

print(output.Data)
print(results)



