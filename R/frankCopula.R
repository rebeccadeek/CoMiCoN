#' The Frank Copula
#'
#' Frank copula function.
#'
#' @param u vector of uniform random variables.
#' @param v vector of uniform random variables.
#' @param theta copula dependence parameter.
#' @param prec the maximal precision to be used in bits.
#'
#' @return \code{frankCopula} returns the frank copula distribution function.
#'
#' @export
#'
#' @examples
#' frankCopula(u = 0.35, v = 0.5, theta = 2)
frankCopula = function(u, v, theta, prec = 200){
  if(class(u)!="numeric" | all(u >= 0) == FALSE | all(u <= 1) == FALSE) {
    stop("ERROR: u must be a numeric vector of Uniform(0,1) random variables.")
  } else if(class(v) !="numeric" | all(v >= 0) == FALSE | all(v <= 1) == FALSE) {
    stop("ERROR: v must be a numeric vector of Uniform(0,1) random variables.")
  } else if(class(theta) !="numeric") {
    stop("ERROR: theta must be a numeric variable.")
  }

  a = Rmpfr::mpfr(exp(-theta), prec) # parameterize a = exp(-theta)

  Cuv = (-1/theta)*as.numeric(log(1 + (((a^u - 1) * (a^v - 1))/(a - 1))))

  return(Cuv)
}
