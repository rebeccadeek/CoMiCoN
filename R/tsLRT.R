#' Two-stage Likelihood Ratio Test
#'
#' Two-stage likelihood ratio test for \eqn{\theta}.
#'
#' @param response data frame with two columns for bivariate pair (x1, x2).
#' @param p1 two-stage estimate of zero-inflation probability for x1.
#' @param mu1 two-stage estimate of beta mean parameter for x1.
#' @param phi1 two-stage estimate of beta dispersion parameter for x1.
#' @param p2 two-stage estimate of zero-inflation probability for x2.
#' @param mu2 two-stage estimate of beta mean parameter for x2.
#' @param phi2 two-stage estimate of beta dispersion parameter for x2.
#' @param theta two-stage estimate of copula dependence parameter.
#' @param thetaVar jackknife variance of \eqn{\theta}.
#' @param thetaNull value of \eqn{\theta} under the null hypothesis. Must be 0 for this implementation of CoMiCoN.
#'
#' @return \code{lrtTheta} returns the two-stage likelihood ratio test statistic.
#'
#' @export
lrtTheta <- function(response, p1, mu1, phi1, p2, mu2, phi2, theta, thetaVar = NULL, thetaNull = 0){
  if(!("data.frame" %in% class(response)) | ncol(response) != 2) {
    stop("ERROR: x must be a data frame with two columns containing the bivariate relative abundances.")
  # } else if(class(lb)!="numeric" | class(ub)!="numeric") {
  #   stop("ERROR: lower and upper must be numeric end points of the interval to be searched.")
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
  } else if(thetaNull != 0) {
    stop("ERROR: This version of CoMiCoN can only handle independence (theta = 0) testing.")
  }

  # calculate the log-likelihood under the null (theta = 0)
  loglik.null <- frankLogLik(theta = thetaNull,
                             p1 = p1, mu1 = mu1, phi1 = phi1,
                             p2 = p2, mu2 = mu2, phi2 = phi2,
                             x = response)

  # calculate the log-likelihood under the alternative (theta = MLE)
  loglik.alt <- frankLogLik(theta = theta,
                            p1 = p1, mu1 = mu1, phi1 = phi1,
                            p2 = p2, mu2 = mu2, phi2 = phi2,
                            x = response)

  # IFF thetaNull != 0: calculate scaling factor
  if(thetaNull != 0){
    # sample size
    n = nrow(response)

    jkV.theta = n*thetaVar

    ## approximate two-stage information for theta
    Idd = (-1/n)*maxLik::hessian(maxLik::maxLik(logLik = frankLogLik,
                                                start = c(theta = theta), x = response,
                                                p1 = p1, mu1 = mu1, phi1 = phi1,
                                                p2 = p2, mu2 = mu2, phi2 = phi2))

    # LRT test statistic
    lr.ts <- -2*(loglik.null - loglik.alt)*(jkV.theta*Idd)^-1
  }

  # else (thetaNull == 0): scaling factor reduces to 1
  else lr.ts <- -2*(loglik.null - loglik.alt)

  return(as.numeric(lr.ts))
}
