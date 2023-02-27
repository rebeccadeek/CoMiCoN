# Functions to find marginal distribution and dependence parameter MLEs

## Marginal MLEs
# Find the regression MLEs, fitted parameters p/mu/phi and shape parameters p/alpha/beta 
# from a zero-inflated beta regression model 
marginalMLE = function(formula.mu, formula.phi = ~ 1, formula.p = ~ 1, df){
  # make inputted formula of class "formula"
  formula.mu = as.formula(formula.mu); formula.phi = as.formula(formula.phi); formula.p = as.formula(formula.p)
  
  # fit ZI-Beta regression model
  model = suppressWarnings( # suppress function warnings
    gamlss::gamlss(formula = formula.mu, sigma.formula = formula.phi, nu.formula = formula.p, 
                   family = gamlss.dist::BEZI, data = df, 
                   control = gamlss::gamlss.control(n.cyc = 40, trace=FALSE))
  )
  
  # check convergence
  convergenceFail = model$converged==FALSE 
  
  # if model failed to converge refit
  if(convergenceFail){
    
    message("Model failed to converge. Refitting.")
    
    model = gamlss::refit(model)
    
    if(model$converged==FALSE){
      message("Model refit and failed to converged.")
    }
  }
  
  # combine estimated regression coefficients 
  zibrCoef = unlist(gamlss::coefAll(model))
  names(zibrCoef) = janitor::make_clean_names(names(zibrCoef))
  
  # estimated P(x=1), mean, and dispersion
  p = fitted(model,"nu"); mu = fitted(model,"mu"); phi = fitted(model,"sigma")
  
  # estimated shape parameters
  alpha = mu*phi; beta = phi - mu*phi
  
  # output
  out = list("coefficients" = zibrCoef, 
             "fitted" = data.frame(p=p, mu=mu, phi=phi), 
             "shape" = data.frame(p=p, alpha=alpha, beta=beta))
  
  return(out)
}

## Theta MLE
# Function to find the MLE of theta using optimizeR function from Rmpfr
# The Rmpfr packages allow for methods of arithmetic functions for arbitrary 
# precision floating point numbers, here often called "mpfr - numbers"
# Cannot find the MLE of theta using optimize() function -- gives error for large theta due to underflow
optTheta <- function(x, lower, upper, p1, alpha1, beta1, p2, alpha2, beta2, prec = 200) {
  
  # define x1 and x2
  x1 <- x[,1]; x2 <- x[,2]
  
  # subset obs. by each combo
  s1 <- which(x1 != 0 & x2 != 0) # number of x1 != 0 and x2 != 0
  s2 <- which(x1 == 0 & x2 != 0) # number of x1 == 0 and x2 != 0
  s3 <- which(x1 != 0 & x2 == 0) # number of x1 != 0 and x2 == 0
  s4 <- which(x1 == 0 & x2 == 0) # number of x1 == 0 and x2 == 0
  
  # source the distribution functions
  source(here::here("R/00b_20200304_distributionFunctions.R"), local = TRUE)
  
  # CDF for each x_i
  F1 <- Fx(x1, p1, alpha1, beta1); F2 <- Fx(x2, p2, alpha2, beta2)
  
  # Function to optimize = log-likelihood; optimize wrt to theta
  f <- function(theta, precision = prec){
    a = Rmpfr::mpfr(exp(-theta), precision) # parameterize a = exp(-theta)
    
    ll <- 
      sum(log((-theta * (a - 1) * a^(F1[s1] + F2[s1])) / ((a^F1[s1] - 1) * (a^F2[s1] - 1) + (a - 1))^2)) + 
      sum(log(((a^F1[s2] - 1) * (a^F2[s2])) / ((a^F1[s2] - 1) * (a^F2[s2] - 1) + (a - 1)))) + 
      sum(log(((a^F2[s3] - 1) * (a^F1[s3])) / ((a^F1[s3] - 1) * (a^F2[s3] - 1) + (a - 1)))) + 
      sum(log((-1/theta) * log(1 + (((a^F1[s4] - 1) * (a^F2[s4] - 1)) / (a - 1)))))
  }
  
  # Numerical analysis to find MLE of theta
  thetaML <- Rmpfr::optimizeR(f, lower = lower, upper = upper, 
                              method = "Brent", maximum = TRUE)$maximum
  
  return(as.numeric(thetaML))
}

## Combine marginal and dependence MLE functions into 1
twoStageMLE = function(formula.mu, formula.phi, formula.p, df, lb = -100, ub = 100, prec = 200){
  
  # marginal MLEs
  marginalMLEs1 = marginalMLE(formula.mu[[1]], formula.phi[[1]], formula.p[[1]], df)
  # rename
  names(marginalMLEs1$coefficients) = gsub("mu", "mu1", names(marginalMLEs1$coefficients))
  names(marginalMLEs1$coefficients) = gsub("sigma", "sigma1", names(marginalMLEs1$coefficients))
  names(marginalMLEs1$coefficients) = gsub("nu", "nu1", names(marginalMLEs1$coefficients))
  names(marginalMLEs1$fitted) = paste0(names(marginalMLEs1$fitted), 1)
  names(marginalMLEs1$shape) = paste0(names(marginalMLEs1$shape), 1)
  
  marginalMLEs2 = marginalMLE(formula.mu[[2]], formula.phi[[2]], formula.p[[2]], df)
  # rename
  names(marginalMLEs2$coefficients) = gsub("mu", "mu2", names(marginalMLEs2$coefficients))
  names(marginalMLEs2$coefficients) = gsub("sigma", "sigma2", names(marginalMLEs2$coefficients))
  names(marginalMLEs2$coefficients) = gsub("nu", "nu2", names(marginalMLEs2$coefficients))
  names(marginalMLEs2$fitted) = paste0(names(marginalMLEs2$fitted), 2)
  names(marginalMLEs2$shape) = paste0(names(marginalMLEs2$shape), 2)
  
  # extract outcome ("x") from each formula.mu
  x = data.frame(x1 = eval(formula.mu[[1]][[2]], df),  # second component of the formula = outcome
                 x2 = eval(formula.mu[[2]][[2]], df) ) # eval(., df) extracts the outcome column from df
  
  thetaMLE = optTheta(x = x, lower = lb, upper = ub, 
                      p1 = marginalMLEs1$shape$p1, alpha1 = marginalMLEs1$shape$alpha1, beta1 = marginalMLEs1$shape$beta1, 
                      p2 = marginalMLEs2$shape$p2, alpha2 = marginalMLEs2$shape$alpha2, beta2 = marginalMLEs2$shape$beta2, 
                      prec = prec)
  
  # merge & output
  MLEs = list(eta = c(marginalMLEs1$coefficients, marginalMLEs2$coefficients, theta=thetaMLE), 
              eta.fitted = cbind(marginalMLEs1$fitted, marginalMLEs2$fitted, theta = thetaMLE), 
              eta.shape = cbind(marginalMLEs1$shape, marginalMLEs2$shape, theta = thetaMLE) )
  
  return(MLEs)
}

# Function to estimate the two-stage covariance matrix using jackknife
# Justification for this approximation in Multivariate Models and Dependence Concepts (Joe) p302
jackknifeVar = function(formula.mu, formula.phi, formula.p, df, eta){
  N = nrow(df) # number of obs.
  
  loo.tsmle = lapply(1:N, function(i) twoStageMLE(formula.mu, formula.phi, formula.p, df[-i,])) # LOO tsMLE
  loo.eta   = lapply(loo.tsmle, function(ls) ls$eta)                                           # extract eta
  loo.op    = lapply(loo.eta, function(eta.i, eta) (eta.i - eta) %o% (eta.i - eta), eta = eta) # outer product
  V = Reduce(`+`, loo.op)
  
  return(V)
}
