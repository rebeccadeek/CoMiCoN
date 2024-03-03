#' Marginal Maximum Likelihood Estimation
#'
#' Maximum likelihood estimation of marginal distribution function parameters.
#'
#' @param formula.mu an object of class formula: a symbolic description of the model to be fitted for the beta mean parameter.
#' @param formula.phi an object of class formula: a symbolic description of the model to be fitted for the beta dispersion parameter.
#' @param formula.p an object of class formula: a symbolic description of the model to be fitted for the zero inflation probability parameter.
#' @param df data frame of outcomes and covariates.
#'
#' @return \code{marginalMLE} returns a list of the following elements:
#' \item{coefficients}{vector of regression coefficients from the models for p, mu, and phi.}
#' \item{fitted}{data frame of fitted values for p, mu, and phi.}
#'
#' @export
marginalMLE = function(formula.mu, formula.phi = ~ 1, formula.p = ~ 1, df){
  if(!("data.frame" %in% class(df))) {
    stop("ERROR: df must be a data frame of covariates for the marginal regression models.")
  }

  # make inputted formula of class "formula"
  formula.mu = stats::as.formula(formula.mu)
  formula.phi = stats::as.formula(formula.phi)
  formula.p = stats::as.formula(formula.p)

  # fit ZI-Beta regression model
  model = suppressWarnings( # suppress function warnings
    gamlss::gamlss(formula = formula.mu, sigma.formula = formula.phi, nu.formula = formula.p,
                   family = gamlss.dist::BEZI, data = df,
                   control = gamlss::gamlss.control(n.cyc = 40, trace=FALSE))
  )

  # check convergence
  convergenceFail = model$converged==FALSE

  # if model failed to converge refit
  if(convergenceFail){

    message("Model failed to converge. Refitting.")

    model = gamlss::refit(model)

    if(model$converged==FALSE){
      message("Model refit and failed to converged.")
    } else {
      message("Model successfully refit.")
    }
  }

  # combine estimated regression coefficients
  zibrCoef = unlist(gamlss::coefAll(model))
  names(zibrCoef) = janitor::make_clean_names(names(zibrCoef))

  # estimated zero probability, mean, and dispersion
  p = stats::fitted(model,"nu"); mu = stats::fitted(model,"mu"); phi = stats::fitted(model,"sigma")

  # output
  out = list("coefficients" = zibrCoef,
             "fitted" = data.frame(p=p, mu=mu, phi=phi)
             )

  return(out)
}
