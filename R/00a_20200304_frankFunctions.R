# Scalar and vectorized functions for:
#   
# 1. Frank copula
# 2. Frank conditional copula (V|U=u and U|V=v)
# 3. Frank copula density
# 
# Note: `mpfr(x, precBits)` from the `Rmprf` package creates multiple (i.e. typically high) 
# precision numbers, to be used in arithmetic and mathematical computations, where x is the 
# high precision number and precBits is maximal precision to be used in bits. `as.numeric()` 
# is used the returned object from mpfr class to numeric.

# Frank Copula 
frankCopula <- function(u, v, theta, prec = 200){
  a = Rmpfr::mpfr(exp(-theta), prec) # parameterize a = exp(-theta)
  
  Cuv = (-1/theta)*as.numeric(log(1 + (((a^u - 1) * (a^v - 1))/(a - 1))))
  
  return(Cuv)
}

# Frank conditional copula, conditional on U
frankConditionalU <- function(u, v, theta, prec = 200){
  a = Rmpfr::mpfr(exp(-theta), prec) # parameterize a = exp(-theta)
  
  C1uv = as.numeric((a^u * (a^v - 1)) / ((a^u - 1) * (a^v - 1) + (a - 1)))
  
  return(C1uv)
}

# Frank conditional copula, conditional on V
frankConditionalV <- function(u, v, theta, prec = 200){
  a = Rmpfr::mpfr(exp(-theta), prec) # parameterize a = exp(-theta)
  
  C1vu = as.numeric((a^v * (a^u - 1)) / ((a^u - 1) * (a^v - 1) + (a - 1)))
  
  return(C1vu)
}

# Frank copula density
frankDensity <- function(u, v, theta, prec = 200){
  a = Rmpfr::mpfr(exp(-theta), prec) # parameterize a = exp(-theta)
  
  cuv = as.numeric((-theta * (a - 1) * a^(u + v)) / ((a^u - 1) * (a^v - 1) + (a - 1))^2)
  
  return(cuv)
}