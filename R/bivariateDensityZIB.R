#' Bivariate Zero-Inflated Beta Density
#'
#' Bivariate frank copula density function with zero-inflated beta margins.
#'
#' @param x1 first vector of relative abundances in the (x1, x2) pair.
#' @param p1 zero-inflation probability for x1.
#' @param mu1 beta mean parameter for x1.
#' @param phi1 beta dispersion parameter for x1.
#' @param x2 second vector of relative abundances in the (x1, x2) pair.
#' @param p2 zero-inflation probability for x2.
#' @param mu2 beta mean parameter for x2.
#' @param phi2 beta dispersion parameter for x2.
#' @param theta copula dependence parameter.
#' @param prec the maximal precision to be used in bits.
#'
#' @return \code{bivariateDensityZIB} returns the bivariate zero-inflated beta density function using a frank copula.
#'
#' @export
bivariateDensityZIB = function(x1, p1, mu1, phi1, x2, p2, mu2, phi2, theta, prec = 200){
  if(class(x1)!="numeric" | all(x1 >= 0) == FALSE | all(x1 < 1) == FALSE |
     class(x2)!="numeric" | all(x2 >= 0) == FALSE | all(x2 < 1) == FALSE) {
    stop("ERROR: x1 and x2 must be numeric vectors of zero-inflated beta random variables with range [0,1).")
  } else if(class(p1) !="numeric" | all(p1 >= 0) == FALSE | all(p1 < 1) == FALSE |
            class(p2) !="numeric" | all(p2 >= 0) == FALSE | all(p2 < 1) == FALSE) {
    stop("ERROR: p1 and p2 must be numeric vectors of zero-inflation probabilities with range [0,1).")
  } else if(class(mu1) !="numeric"  | all(mu1 >= 0) == FALSE | all(mu1 < 1) == FALSE |
            class(mu2) !="numeric"  | all(mu2 >= 0) == FALSE | all(mu2 < 1) == FALSE) {
    stop("ERROR: mu1 and mu2 must be numeric vectors of mean parameters with range [0,1).")
  } else if(class(phi1) !="numeric" | all(phi1 > 0) == FALSE |
            class(phi2) !="numeric" | all(phi2 > 0) == FALSE) {
    stop("ERROR: phi1 and phi2 must be positive and non-zero numeric vectors of dispersion parameters.")
  } else if(class(theta) !="numeric") {
    stop("ERROR: theta must be a numeric variable.")
  }

  # pdf
  fx1 = dzib(x1, p1, mu1, phi1) # density of x1
  fx2 = dzib(x2, p2, mu2, phi2) # density of x2

  # cdf
  Fx1 = pzib(x1, p1, mu1, phi1) # cdf at x1
  Fx2 = pzib(x2, p2, mu2, phi2) # cdf at x2

  # empty joint distribution vector
  fx1x2 = vector(length = length(x1))

  # subset obs. by each combo
  s1 = which(x1 != 0 & x2 != 0) # number of x1 != 0 & x2 != 0
  s2 = which(x1 == 0 & x2 != 0) # number of x1 == 0 & x2 != 0
  s3 = which(x1 != 0 & x2 == 0) # number of x1 != 0 & x2 == 0
  s4 = which(x1 == 0 & x2 == 0) # number of x1 == 0 & x2 == 0

  # joint distribution -- using frank copula
  fx1x2[s1] = frankDensity(Fx1[s1], Fx2[s1], theta)*fx1[s1]*fx2[s1] # joint density for s1
  fx1x2[s2] = fx2[s2]*frankConditionalV(Fx1[s2], Fx2[s2], theta)    # joint density for s2
  fx1x2[s3] = fx1[s3]*frankConditionalU(Fx1[s3], Fx2[s3], theta)    # joint density for s3
  fx1x2[s4] = frankCopula(Fx1[s4], Fx2[s4], theta)                  # joint density for s4

  return(fx1x2)
}
