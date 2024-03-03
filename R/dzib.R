#' The Zero-Inflated Beta Distribution
#'
#' Density function for the zero-inflated beta distribution with zero-inflation probability equal to \code{p}, mean equal to \code{mu} and dispersion equal to \code{phi}.
#'
#' @param x vector of quantiles.
#' @param p vector of zero-inflation probabilities.
#' @param mu vector of means.
#' @param phi vector of dispersions.
#'
#' @return \code{dzib} returns the density.
#'
#' @export
#'
#' @examples
#' dzib(x = 0.8, p = 0.2, mu = 0.5, phi = 0.1)
dzib = function(x, p, mu, phi){
  if(class(x)!="numeric" | all(x >= 0) == FALSE | all(x < 1) == FALSE) {
    stop("ERROR: x must be a numeric vector of zero-inflated beta random variables with range [0,1).")
  } else if(class(p) !="numeric" | all(p >= 0) == FALSE | all(p < 1) == FALSE) {
    stop("ERROR: p must be a numeric vector of zero-inflation probabilities with range [0,1).")
  } else if(class(mu) !="numeric"  | all(mu >= 0) == FALSE | all(mu < 1) == FALSE) {
    stop("ERROR: mu must be a numeric vector of mean parameters with range [0,1).")
  } else if(class(phi) !="numeric" | all(phi > 0) == FALSE) {
    stop("ERROR: phi must be a positive and non-zero numeric vector of dispersion parameters.")
  }

  # from mean/dispersion to shape parameters
  alpha = mu*phi; beta = phi - mu*phi

  fx = ifelse(x == 0, p, (1 - p) * stats::dbeta(x, alpha, beta))

  return(fx)
}
