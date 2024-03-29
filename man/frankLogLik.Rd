% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/frankLogLik.R
\name{frankLogLik}
\alias{frankLogLik}
\title{Log-likelihood for the Bivariate Zero-Inflated Beta Density}
\usage{
frankLogLik(theta, p1, mu1, phi1, p2, mu2, phi2, x, prec = 200)
}
\arguments{
\item{theta}{copula dependence parameter.}

\item{p1}{zero-inflation probability for x1.}

\item{mu1}{beta mean parameter for x1.}

\item{phi1}{beta dispersion parameter for x1.}

\item{p2}{zero-inflation probability for x2.}

\item{mu2}{beta mean parameter for x2.}

\item{phi2}{beta dispersion parameter for x2.}

\item{x}{data frame with two columns for bivariate pair (x1, x2).}

\item{prec}{the maximal precision to be used in bits.}
}
\value{
Two-stage estimate of \eqn{\theta}.
}
\description{
Log-likelihood for the bivariate frank copula density function with a zero-inflated beta margins.
}
