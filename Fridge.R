#' Fridge
#'
#' Function for calculating the fridge-ridge tuning parameter and prediction, with plot of MSE curve  
#'
#' @param X Matrix with covariates of dimension n x p
#' @param y Vector of outcomes of length n
#' @param x0 Matrix of focus covariates of length p x t. Each column is one
#' focus covariate.
#' @param plug.in Character specifying the plug.in estimate, possibly 'OLS' (ordinary least squares) and 'RLO0CV' (ridge with leave-one-out cross-validation). 
#' The OLS option can only be used for p smaller than n. 
#' @param plot.curve Logical. Should the MSE curve and minimum value be ploted? 
#' @details The minimization routine is carried out with several different start values ranging from 1 to 10^7. If plot.curve = TRUE, the global minimum
#' is ploted. 
#' @return Result, a list with the following 
#' \item{fridge.tuning}{Numerical, the fridge tuning parameter} 
#' \item{loocv.tuning}{Numerical, the leave-one-out cross-validation tuning parameter}
#' \item{fridge.pred}{Numerical, the fridge prediction} 
#' \item{loocv.pred}{Numerical, the leave-one-out cross-validation prediction}
#' @keywords Focused tuning
#' @export
#' @examples
#' fridge()

fridge <- function(X,y,x0,plug.in = c('OLS','RLOOCV'),plot.curve=TRUE){
  
  if(dim(X)[1]!=length(y)){stop('Dimensions of y and X do not agree');}
  if(dim(X)[2]!=dim(x0)[1]){stop('Dimensions of x0 and X do not agree');}
  if(dim(X)[1]<dim(X)[2]&&plug.in=='OLS'){stop('OLS plug-in cannot be used for high-dimensional data');}
  
  result <- list()
  p <- dim(X)[2]
  n <- dim(X)[1]
  num.focus.cov <- dim(x0)[2]
  
  ### Plug-in estimator (indepenedent of focus covariate)
  
  start.time <- Sys.time()  # start time measurement
  
  svd <- svd(X)
  opt.loocv <- optim(par = 1,tuning.loocv,y=y,svd=svd,n=n,lower=1e-15,upper = Inf, method = "L-BFGS-B",control = list(factr=1e4))$par
  
  #Decide on a set of start values
  start <- 10^seq(1,7,1)
  value <- c(); param <- c();
  
  # Calculate plug-in estimates
  if(plug.in=='OLS'){
    sigma2_hat <- sum((y - svd$u%*%t(svd$u)%*%y)^2)/(n-p)
  } else {
    sigma2_hat <- sum((y - svd$u%*%diag(svd$d^2/(svd$d^2 + opt.loocv))%*%t(svd$u)%*%y)^2)/(n-sum(diag(svd$u%*%diag(svd$d^2/(svd$d^2 + opt.loocv))%*%t(svd$u))))
  }  
  
  # Compute loocv tuning parameter and prediction
  result$loocv.tuning <- opt.loocv
  result$loocv.beta <- svd$v%*%diag(svd$d/(svd$d^2  + opt.loocv))%*%t(svd$u)%*% y
  result$loocv.pred <- t(x0)%*%result$loocv.beta
  
  # run time measurement
  time.plug.in <- as.numeric(difftime(Sys.time(), start.time, units = "sec"))
  
  
  ### Tuning parameter calibration for focus covariate
  
  # run time measurement of the focus covariant dependent computation
  start.time <- Sys.time()
  
  fridge.tuning <- rep(NA, num.focus.cov)
  fridge.beta <- matrix(NA, nrow=p, ncol=num.focus.cov)
  fridge.pred <- rep(NA, num.focus.cov)
  
  for(t in 1:num.focus.cov) {
    z <- x0[, t]  # focus covariate
    for(l in 1:length(start)) {
      if(plug.in=='OLS'){
        optimfridge <- optim(par = start[l],tuning.fridge.ols.sim,y=y,svd = svd,sigma2_hat=sigma2_hat,x0=z,lower=0,upper = Inf, method = "L-BFGS-B",control = list(factr=1e2))
      } else {
        optimfridge <- optim(par = start[l],tuning.fridge.ridge.sim,y=y,svd = svd,sigma2_hat=sigma2_hat,cv.tuning=opt.loocv,x0=z,lower=0,upper = Inf, method = "L-BFGS-B",control = list(factr=1e2))
      } 
      value[l] <- optimfridge$value   
      param[l] <- optimfridge$par  
    }
    opt.fic  <- param[which.min(value)] 
    
    # if(plot.curve==TRUE){
    #   end.point <- opt.fic*2*(opt.fic<10^5) + 10^5*(opt.fic>10^5) 
    #   if(plug.in=='OLS'){
    #     curve(sapply(x,tuning.fridge.ols.sim,y=y,svd = svd,sigma2_hat=sigma2_hat,x0=z),0,end.point,n=500,
    #           main= paste0('The minimum MSE is given at ',round(opt.fic,1)),xlab='Tuning parameter',ylab='MSE')
    #     abline(v=opt.fic,col=2)
    #   } else {
    #     curve(sapply(x,tuning.fridge.ridge.sim,y=y,svd = svd,sigma2_hat=sigma2_hat,cv.tuning=opt.loocv,x0=z),0,end.point,n=500,
    #           main= paste0('The minimum MSE is given at ',round(opt.fic,1)),xlab='Tuning parameter',ylab='MSE')
    #     abline(v=opt.fic,col=2)
    #   }  
    # }  
    
    # optimal tuning parameter
    fridge.tuning[t] <- opt.fic
    
    # ridge prediction
    fridge.beta[, t] <- svd$v%*%diag(svd$d/(svd$d^2  + opt.fic))%*%t(svd$u)%*% y
    fridge.pred[t] <- t(z)%*%fridge.beta[, t]
  }
  
  # run time measurement
  time.focus.cov <- as.numeric(difftime(Sys.time(), start.time, units = "sec"))
  time.focus.cov <- time.focus.cov / num.focus.cov  # mean run time
  
  # save and return results
  result$fridge.tuning <- fridge.tuning
  result$fridge.pred <- fridge.pred
  result$fridge.beta <- fridge.beta
  result$time.plug.in <- time.plug.in
  result$time.focus.cov <- time.focus.cov
  return(result)
}

#' Cross-validation
#'
#' Function for the ridge regression leave-one-out cross-validation error
#'
#' @param svd SVD of the n x p data matrix
#' @param n Number of observations
#' @param lambda Ridge tuning parameter 
#' @details The function gives the explitic expression of the LOOCV error for ridge regression
#' @keywords Focused tuning
#' @export
#' @examples
#' tuning.loocv

tuning.loocv <- function(lambda,svd,n,y){
  mean((((diag(n)-svd$u %*% diag(svd$d^2/(svd$d^2 +lambda)) %*% t(svd$u)) %*% y)/(1-diag(svd$u %*% diag(svd$d^2/(svd$d^2 +lambda)) %*% t(svd$u))))^2)
}

#' Simplified fridge-ridge 
#'
#' Function for fridge-ridge error 
#'
#' @param svd SVD of the n x p data matrix
#' @param lambda Ridge tuning parameter 
#' @param x0 Focus covariate
#' @param sigma2_hat Plug-in estimate of the error variance
#' @param init.tuning Initial ridge tuning parameter
#' @details The function gives error for fridge-ridge estimator
#' @keywords Focused tuning
#' @export
#' @examples
#' tuning.fridge.ridge.sim

tuning.fridge.ridge.sim <- function(lambda,svd,y,sigma2_hat,x0,cv.tuning){
  return((lambda*t(x0)%*%svd$v %*% diag(1/(svd$d^2 + lambda)*(svd$d)/(svd$d^2+cv.tuning)) %*% t(svd$u) %*% y)^2
         +  sigma2_hat*t(x0)%*%svd$v%*% diag(svd$d^2/(svd$d^2 + lambda)^2) %*% t(svd$v)%*%x0)
}

#' Simplified fridge-OLS 
#'
#' Function for fridge-OLS error 
#'
#' @param svd SVD of the n x p data matrix
#' @param lambda Ridge tuning parameter 
#' @param x0 Focus covariate
#' @param sigma2_hat Plug-in estimate of the error variance
#' @details Can only be used in the low-dimensional case, for n < p
#' @keywords Focused tuning
#' @export
#' @examples
#' tuning.fridge.ols.sim

tuning.fridge.ols.sim <- function(lambda,svd,sigma2_hat,x0,y){
  return((lambda*t(x0)%*%svd$v %*% diag(1/(svd$d^2 + lambda)*1/svd$d) %*% t(svd$u) %*% y)^2
         +  sigma2_hat*t(x0)%*%svd$v%*% diag(svd$d^2/(svd$d^2 + lambda)^2) %*% t(svd$v)%*%x0)
}




