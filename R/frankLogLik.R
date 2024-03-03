#' Log-likelihood for the Bivariate Zero-Inflated Beta Density
#'
#' Log-likelihood for the bivariate frank copula density function with a zero-inflated beta margins.
#'
#' @param theta copula dependence parameter.
#' @param p1 zero-inflation probability for x1.
#' @param mu1 beta mean parameter for x1.
#' @param phi1 beta dispersion parameter for x1.
#' @param p2 zero-inflation probability for x2.
#' @param mu2 beta mean parameter for x2.
#' @param phi2 beta dispersion parameter for x2.
#' @param x data frame with two columns for bivariate pair (x1, x2).
#' @param prec the maximal precision to be used in bits.
#'
#' @return Two-stage estimate of \eqn{\theta}.
#'
#' @export
frankLogLik <- function(theta, p1, mu1, phi1, p2, mu2, phi2, x, prec = 200) {
  if(class(theta) !="numeric") {
    stop("ERROR: theta must be a numeric variable.")
  } else if(class(p1) !="numeric" | all(p1 >= 0) == FALSE | all(p1 < 1) == FALSE |
            class(p2) !="numeric" | all(p2 >= 0) == FALSE | all(p2 < 1) == FALSE) {
    stop("ERROR: p1 and p2 must be numeric vectors of zero-inflation probabilities with range [0,1).")
  } else if(class(mu1) !="numeric"  | all(mu1 >= 0) == FALSE | all(mu1 < 1) == FALSE |
            class(mu2) !="numeric"  | all(mu2 >= 0) == FALSE | all(mu2 < 1) == FALSE) {
    stop("ERROR: mu1 and mu2 must be numeric vectors of mean parameters with range [0,1).")
  } else if(class(phi1) !="numeric" | all(phi1 > 0) == FALSE |
            class(phi2) !="numeric" | all(phi2 > 0) == FALSE) {
    stop("ERROR: phi1 and phi2 must be positive and non-zero numeric vectors of dispersion parameters.")
  } else if(!("data.frame" %in% class(x)) | ncol(x) != 2) {
    stop("ERROR: x must be a data frame with two columns containing the bivariate relative abundances.")
  }

  # parameterize a = exp(-theta)
  a = Rmpfr::mpfr(exp(-theta), prec)

  # define x1 and x2
  x1 <- x[,1]; x2 <- x[,2]

  # Check if theta != 0
  thetaCheck <- theta != 0

  if(thetaCheck){ # theta != 0, use copula log-likelihood
    ll = sum(log(bivariateDensityZIB(x1, p1, mu1, phi1, x2, p2, mu2, phi2, theta, prec = 200)))
  }

  else { # theta == 0, model reduce to independence-- product of marginal pdfs
    # PDF for each x_i
    f1 = dzib(x1, p1, mu1, phi1); f2 = dzib(x2, p2, mu2, phi2)

    ll = sum(log(f1*f2))
  }

  return(as.numeric(ll))
}
