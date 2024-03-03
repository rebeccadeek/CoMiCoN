#' Two-stage Estimation of \eqn{\theta}
#'
#' Two-stage estimation of copula dependence parameter \eqn{\theta}.
#'
#' @param x data frame with two columns for bivariate pair (x1, x2).
#' @param lower the lower end point of the interval to be searched in univariate optimizer.
#' @param upper the lower end point of the interval to be searched in univariate optimizer.
#' @param p1 two-stage estimate of the zero-inflation probability for x1.
#' @param mu1 two-stage estimate of the beta mean parameter for x1.
#' @param phi1 two-stage estimate of the beta dispersion parameter for x1.
#' @param p2 two-stage estimate of the zero-inflation probability for x2.
#' @param mu2 two-stage estimate of the beta mean parameter for x2.
#' @param phi2 two-stage estimate of the beta dispersion parameter for x2.
#' @param prec the maximal precision to be used in bits.
#'
#' @return \code{thetaTSMLE} returns the two-stage estimate of \eqn{\theta}.
#'
#' @export
thetaTSMLE <- function(x, lower, upper, p1, mu1, phi1, p2, mu2, phi2, prec = 200) {
  if(!("data.frame" %in% class(x)) | ncol(x) != 2) {
    stop("ERROR: x must be a data frame with two columns containing the bivariate relative abundances.")
  } else if(class(lower)!="numeric" | class(upper)!="numeric") {
    stop("ERROR: lower and upper must be numeric end points of the interval to be searched.")
  } else if(class(p1) !="numeric" | all(p1 >= 0) == FALSE | all(p1 < 1) == FALSE |
            class(p2) !="numeric" | all(p2 >= 0) == FALSE | all(p2 < 1) == FALSE) {
    stop("ERROR: p1 and p2 must be numeric vectors of zero-inflation probabilities with range [0,1).")
  } else if(class(mu1) !="numeric"  | all(mu1 >= 0) == FALSE | all(mu1 < 1) == FALSE |
            class(mu2) !="numeric"  | all(mu2 >= 0) == FALSE | all(mu2 < 1) == FALSE) {
    stop("ERROR: mu1 and mu2 must be numeric vectors of mean parameters with range [0,1).")
  } else if(class(phi1) !="numeric" | all(phi1 > 0) == FALSE |
            class(phi2) !="numeric" | all(phi2 > 0) == FALSE) {
    stop("ERROR: phi1 and phi2 must be positive and non-zero numeric vectors of dispersion parameters.")
  }

  # define x1 and x2
  x1 <- x[,1]; x2 <- x[,2]

  # subset obs. by each combo
  s1 <- which(x1 != 0 & x2 != 0) # number of x1 != 0 and x2 != 0
  s2 <- which(x1 == 0 & x2 != 0) # number of x1 == 0 and x2 != 0
  s3 <- which(x1 != 0 & x2 == 0) # number of x1 != 0 and x2 == 0
  s4 <- which(x1 == 0 & x2 == 0) # number of x1 == 0 and x2 == 0

  # CDF for each x_i
  F1 <- pzib(x1, p1, mu1, phi1); F2 <- pzib(x2, p2, mu2, phi2)

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
  thetaTSML <- Rmpfr::optimizeR(f, lower = lower, upper = upper, method = "Brent", maximum = TRUE)$maximum

  return(as.numeric(thetaTSML))
}
