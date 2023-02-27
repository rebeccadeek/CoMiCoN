# Likelihood ratio test for theta ($H_{0}: \theta = 0$ vs. $H_{1}: \theta \neq 0$)
#   1. Simulate data under the null (theta = 0)
#   2. Calculate the MLE of the simulated data
#   3. Calculate the log-likelihood under the MLE value
#   4. Calculate the log-likelihood under the true null (theta = 0)
#   5. Calculate sc = scaling factor for LR test statistic
#   6. -2log(LR)*sc should follow a chi-square 1 distribution under the null
lrtTheta <- function(response, eta.shape, thetaVar, thetaNull = 0, lb = -100, ub = 100){
  # source log-likelihood
  source(here::here("code/01c_20200312_loglikFunction.R"), local = TRUE)
  
  # eta = vector of parameters
  p1 = eta.shape[["p1"]]; alpha1 = eta.shape[["alpha1"]]; beta1 = eta.shape[["beta1"]]
  p2 = eta.shape[["p2"]]; alpha2 = eta.shape[["alpha2"]]; beta2 = eta.shape[["beta2"]]
  theta = eta.shape[["theta"]][[1]]
  
  # calculate the log-likelihood under the null (theta = 0)
  loglik.null <- llFunction(theta = thetaNull, p1 = p1, alpha1 = alpha1, beta1 = beta1, 
                            p2 = p2, alpha2 = alpha2, beta2 = beta2, x = response)
  
  # calculate the log-likelihood under the alternative (theta = MLE)
  loglik.alt <- llFunction(theta = theta, p1 = p1, alpha1 = alpha1, beta1 = beta1, 
                           p2 = p2, alpha2 = alpha2, beta2 = beta2, x = response)
  
  # IFF thetaNull != 0: calculate scaling factor
  if(thetaNull != 0){
    # sample size
    n = nrow(response)
    
    # scaling factor = (v*Idd)^-1
    ## v = n*varTheta
    
    jkV.theta = n*thetaVar
    
    ## approximate two-stage information for theta
    Idd = (-1/n)*maxLik::hessian(maxLik::maxLik(logLik = llFunction, #grad = gt, 
                                                start = c(theta = theta), x = response, 
                                                p1 = p1, alpha1 = alpha1, beta1 = beta1, 
                                                p2 = p2, alpha2 = alpha2, beta2 = beta2))
    
    # LRT test statistic
    lr.ts <- -2*(loglik.null - loglik.alt)*(jkV.theta*Idd)^-1
  }
  
  # else (thetaNull == 0): scaling factor reduces to 1
  else lr.ts <- -2*(loglik.null - loglik.alt)
  
  return(as.numeric(lr.ts))
}