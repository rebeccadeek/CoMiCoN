# Log-likelihood function in terms of the Frank copula
llFunction <- function(theta, p1, alpha1, beta1, p2, alpha2, beta2, x, prec = 200) {
  a = Rmpfr::mpfr(exp(-theta), prec)   # parameterize a = exp(-theta)
  
  # source the distribution functions
  source(here::here("R/00b_20200304_distributionFunctions.R"), local = TRUE)
  
  # define x1 and x2
  x1 <- x[,1]; x2 <- x[,2]
  
  # Check if theta != 0
  thetaCheck <- theta != 0
  
  if(thetaCheck){ # theta != 0, use copula log-likelihood
    ll <- sum(log(fxy(x1, p1, alpha1, beta1, x2, p2, alpha2, beta2, theta, prec = 200)))
  }
  
  else { # theta == 0, model reduce to independence-- product of marginal pdfs
    # PDF for each x_i
    f1 <- fx(x1, p1, alpha1, beta1); f2 <- fx(x2, p2, alpha2, beta2)
    
    ll = sum(log(f1*f2))
  }
  
  return(as.numeric(ll))
}

# Vectorize llFunction wrt theta
llFunctionVect <- Vectorize(llFunction, vectorize.args = "theta")

# Gradient wrt theta
gt = function(theta, p1, alpha1, beta1, p2, alpha2, beta2, x, prec = 200){
  a = Rmpfr::mpfr(exp(-theta), prec)   # parameterize a = exp(-theta)
  
  # define x1 and x2
  x1 <- x[,1]; x2 <- x[,2]
  
  # subset obs. by each combo
  s1 <- which(x1 != 0 & x2 != 0) # number of x1 != 0 and x2 != 0
  s2 <- which(x1 == 0 & x2 != 0) # number of x1 == 0 and x2 != 0
  s3 <- which(x1 != 0 & x2 == 0) # number of x1 != 0 and x2 == 0
  s4 <- which(x1 == 0 & x2 == 0) # number of x1 == 0 and x2 == 0
  
  # number of obs. in each combo
  n1 <- length(s1) # number of x1 != 0 and x2 != 0
  n2 <- length(s2) # number of x1 == 0 and x2 != 0
  n3 <- length(s3) # number of x1 != 0 and x2 == 0
  n4 <- length(s4) # number of x1 == 0 and x2 == 0
  
  # source the distribution functions
  source(here::here("R/00b_20200304_distributionFunctions.R"), local = TRUE)
  
  # CDF for each x_i
  F1 <- Fx(x1, p1, alpha1, beta1); F2 <- Fx(x2, p2, alpha2, beta2)
  
  dlldT <- 
    # likelihood contribution from s1
    n1/theta - 
    n1*a/(a-1) - 
    sum(F1[s1] + F2[s1]) +
    2*sum((F1[s1]*a^F1[s1]*(a^F2[s1] - 1) + F2[s1]*a^F2[s1]*(a^F1[s1] - 1) + a) / ((a^F1[s1] - 1)*(a^F2[s1] - 1) + (a - 1))) - 
    # likelihood contribution from s2
    sum((F1[s2]*a^F1[s2])/(a^F1[s2] - 1)) - 
    sum((F1[s2]*a^F1[s2]*(a^F2[s2] - 1) + F2[s2]*a^F2[s2]*(a^F1[s2] - 1) + a) / ((a^F1[s2] - 1)*(a^F2[s2] - 1) + (a - 1))) - 
    sum(F2[s2]) - 
    # likelihood contribution from s3
    sum((F2[s3]*a^F2[s3])/(a^F2[s3] - 1)) - 
    sum((F1[s3]*a^F1[s3]*(a^F2[s3] - 1) + F2[s3]*a^F2[s3]*(a^F1[s3] - 1) + a) / ((a^F1[s3] - 1)*(a^F2[s3] - 1) + (a - 1))) - 
    sum(F1[s3]) - 
    # likelihood contribution from s4
    n4/theta - 
    sum(((F1[s4]*a^F1[s4]*(a^F2[s4] - 1) - F2[s4]*a^F2[s4]^(a^F1[s4] - 1))*(a - 1) + (a^F1[s4] - 1)*(a^F2[s4] - 1)*a) / 
          (log(1 + ((a^F1[s4] - 1)*(a^F2[s4] - 1))/(a - 1))*(1 + ((a^F1[s4] - 1)*(a^F2[s4] - 1))/(a - 1))*(a - 1)^2))
  
  return(as.numeric(dlldT))
}

# 2nd derivative of ll wrt theta
d2lldTheta = function(theta, p1, alpha1, beta1, p2, alpha2, beta2, x, prec = 200){
  a = Rmpfr::mpfr(exp(-theta), prec)   # parameterize a = exp(-theta)
  eno = Rmpfr::mpfr(exp(-1), prec)     # parameterize eno = exp(-1)
  
  # define x1 and x2
  x1 <- x[,1]; x2 <- x[,2]
  
  # subset obs. by each combo
  s1 <- which(x1 != 0 & x2 != 0) # which x1 != 0 and x2 != 0
  s2 <- which(x1 == 0 & x2 != 0) # which x1 == 0 and x2 != 0
  s3 <- which(x1 != 0 & x2 == 0) # which x1 != 0 and x2 == 0
  s4 <- which(x1 == 0 & x2 == 0) # which x1 == 0 and x2 == 0
  s23 <- sort(c(s2, s3))         # which in s2 or s3
  
  # number of obs. in each combo
  n1 <- length(s1) # number of x1 != 0 and x2 != 0
  n2 <- length(s2) # number of x1 == 0 and x2 != 0
  n3 <- length(s3) # number of x1 != 0 and x2 == 0
  n4 <- length(s4) # number of x1 == 0 and x2 == 0
  
  # source the distribution functions
  source(here::here("R/00b_20200304_distributionFunctions.R"), local = TRUE)
  
  # CDF for each x_i
  F1 <- Fx(x1, p1, alpha1, beta1); F2 <- Fx(x2, p2, alpha2, beta2)
  
  # part of contribution from s4
  B    = 1 + ( ( a^(F1[s4] + F2[s4]) - a^F1[s4] - a^F2[s4] + 1 ) / (a - 1) )
  dBdT = ( ( -(F1[s4] + F2[s4]) * eno^(F1[s4] + F2[s4]) + eno^F1[s4] + eno^F2[s4] ) 
           * (a^-2 - a) + (a^(F1[s4] + F2[s4]) - a^F1[s4] - a^F2[s4] + 1) * a ) / (a - 1)^2
  
  d2ll.dT <- 
    # contribution from s1
    - n1/theta^2 - 
    n1*a/(a-1)^2 - 
    2*sum( ( (F1[s1] + F2[s1])^2 * eno^(F1[s1] + F2[s1]) - F1[s1]^2 * eno^F1[s1] - F2[s1]^2 * eno^F2[s1] + 1 ) / ( eno^(F1[s1] + F2[s1]) - eno^F1[s1] - eno^F2[s1] + 1 ) ) + 
    2*sum( ( ( -(F1[s1] + F2[s1]) * eno^(F1[s1] + F2[s1]) - F1[s1] * eno^F1[s1] - F2[s1] * eno^F2[s1] - 1 ) / ( eno^(F1[s1] + F2[s1]) - eno^F1[s1] - eno^F2[s1] + 1 ) )^2 ) + 
    # contribution from s2
    sum( -(F1[s2])^2 * a^F1[s2] / (a^F1[s2] - 1)^2 ) + 
    # contribution from s3
    sum( -(F2[s3])^2 * a^F2[s3] / (a^F2[s3] - 1)^2 ) + 
    # contribution from s2 and s3
    sum( ( (F1[s23] + F2[s23])^2 * eno^(F1[s23] + F2[s23]) - F1[s23]^2 * eno^F1[s23] - F2[s23]^2 * eno^F2[s23] + 1 ) / ( eno^(F1[s23] + F2[s23]) - eno^F1[s23] - eno^F2[s23] + 1 ) ) - 
    sum( ( ( -(F1[s23] + F2[s23]) * eno^(F1[s23] + F2[s23]) + F1[s23] * eno^F1[s23] + F2[s23] * eno^F2[s23] - 1 ) / ( eno^(F1[s23] + F2[s23]) - eno^F1[s23] - eno^F2[s23] + 1 ) )^2 ) + 
    # contribution from s4
    n4/theta^2 + 
    sum( ( ( -(F1[s4] + F2[s4]) * eno^(F1[s4] + F2[s4]) + F1[s4] * eno^F1[s4] + F2[s4] * eno^F2[s4] ) * (-2*a^2 - a^(-1)) - (-(F1[s4] + F2[s4] + 1) * a^(F1[s4] + F2[s4] + 1) + (F1[s4] + 1) * a^(F1[s4] + 1) + (F2[s4] + 1) * a^(F2[s4] + 1) - a) ) / ( log(B) * B * (a - 1)^2 ) ) - 
    sum( ( ( -(F1[s4] + F2[s4]) * eno^(F1[s4] + F2[s4]) + F1[s4] * eno^F1[s4] + F2[s4] * eno^F2[s4] ) * (-2*a^2 - a) - (a^(F1[s4] + F2[s4] + 1) - a^(F1[s4] + 1) - a^(F2[s4] + 1) + a) ) * ( dBdT * (a - 1)^2 + log(B) * dBdT * (a - 1)^2 - 2 * log(B) * B * (a^2 - a) ) / ( log(B) * B * (a - 1)^2 )^2 )
  
  return(as.numeric(d2ll.dT))
}
