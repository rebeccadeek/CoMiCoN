# Zero-inflated beta marginal pdf and cdf functions
## univariate marginal density of x
fx <- function(x, p, alpha, beta){
  fxTemp = ifelse(x == 0, p, (1 - p) * dbeta(x, alpha, beta))
  
  return(fxTemp)
}

## univariate marginal distribution of x
Fx <- function(x, p, alpha, beta){
  FxTemp = ifelse(x == 0, p, p + (1 - p) * pbeta(x, alpha, beta))
  
  return(FxTemp)
}

# Joint distribution function for f(x1, x2)
fxy <- function(x1, p1, alpha1, beta1, x2, p2, alpha2, beta2, theta, prec = 200){
  # pdf
  fx1 <- fx(x1, p1, alpha1, beta1); fx2 <- fx(x2, p2, alpha2, beta2) # density of x1 and x2
  
  # cdf
  Fx1 <- Fx(x1, p1, alpha1, beta1); Fx2 <- Fx(x2, p2, alpha2, beta2) # cdf at x1 and x2
  
  # source frank copula functions
  source(here::here("R/00a_20200304_frankFunctions.R"), local = TRUE)
  
  # empty joint distribution vector
  fx1x2Temp <- vector(length = length(x1))
  
  # subset obs. by each combo
  s1 <- which(x1 != 0 & x2 != 0) # number of x1 != 0 and x2 != 0
  s2 <- which(x1 == 0 & x2 != 0) # number of x1 == 0 and x2 != 0
  s3 <- which(x1 != 0 & x2 == 0) # number of x1 != 0 and x2 == 0
  s4 <- which(x1 == 0 & x2 == 0) # number of x1 == 0 and x2 == 0
  
  # joint distribution
  fx1x2Temp[s1] <- frankDensity(Fx1[s1], Fx2[s1], theta)*fx1[s1]*fx2[s1] # joint density if x1 != 0 and x2 != 0
  fx1x2Temp[s2] <- fx2[s2]*frankConditionalV(Fx1[s2], Fx2[s2], theta)    # joint density if x1 == 0 and x2 != 0
  fx1x2Temp[s3] <- fx1[s3]*frankConditionalU(Fx1[s3], Fx2[s3], theta)    # joint density if x1 != 0 and x2 == 0
  fx1x2Temp[s4] <- frankCopula(Fx1[s4], Fx2[s4], theta)                  # joint density if x1 == 0 and x2 == 0
  
  return(fx1x2Temp)
}

# Conditional distribution function for f(x2|x1)
# f(x2|x1)
fygx <- function(x1, p1, alpha1, beta1, x2, p2, alpha2, beta2, theta, prec = 200){
  # pmf
  fx1 <- fx(x1, p1, alpha1, beta1); fx2 <- fx(x2, p2, alpha2, beta2) # density of x1 and x2
  
  # cdf
  Fx1 <- Fx(x1, p1, alpha1, beta1); Fx2 <- Fx(x2, p2, alpha2, beta2) # cdf at x1 and x2
  
  # source frank copula functions
  source(here::here("R/00a_20200304_frankFunctions.R"), local = TRUE)
  
  # empty conditional distribution vector
  fx2gx1Temp <- vector(length = length(x1))
  
  # subset obs. by each combo
  s1 <- which(x1 != 0 & x2 != 0) # number of x1 != 0 and x2 != 0
  s2 <- which(x1 == 0 & x2 != 0) # number of x1 == 0 and x2 != 0
  s3 <- which(x1 != 0 & x2 == 0) # number of x1 != 0 and x2 == 0
  s4 <- which(x1 == 0 & x2 == 0) # number of x1 == 0 and x2 == 0
  
  # conditional distribution
  fx2gx1Temp[s1] <- frankDensityVect(Fx1[s1], Fx2[s1], theta)*fx2[s1]                # conditional density if x1 != 0 and x2 != 0
  fx2gx1Temp[s2] <- (fx2[s2]*frankConditionalVVect(Fx1[s2], Fx2[s2], theta))/fx1[s2] # conditional density if x1 == 0 and x2 != 0
  fx2gx1Temp[s3] <- frankConditionalUVect(Fx1[s3], Fx2[s3], theta)                   # conditional density if x1 != 0 and x2 == 0
  fx2gx1Temp[s4] <- frankCopulaVect(Fx1[s4], Fx2[s4], theta)/fx1[s4]                 # conditional density if x1 == 0 and x2 == 0
  
  return(fx2gx1Temp)
}
