#' CoMiCoN
#'
#'  Implementation of Copulas with Mixture Margins for Covariation Networks using two-stage estimation and testing.
#'
#' @param abd data frame of microbial relative abundances. Rows are samples and columns are taxa.
#' @param covars data frame of covariates to adjust for in the marginal regression models. Rows should be in the same order as \code{abd} and rownames should match.
#' @param test logical indicating if two-stage likelihood ratio testing should be performed. Default is \code{FALSE}.
#' @param ncores number of cores for parallelization. Default value is 1 (no parallel computing).
#'
#' @return \code{comicon} returns a data frame with two-stage estimated marginal and dependence parameters (and optional two-stage likelihood ratio test statistic) for all possible pairs of microbes.
#'
#' @export
comicon = function(abd, covars = NULL, test = FALSE, ncores = 1) {
  if(!("data.frame" %in% class(abd))) {
    stop("ERROR: abd must be a data frame of relative abundances.")
  } else if(!is.null(covars) & !("data.frame" %in% class(covars))) {
    stop("ERROR: covars must be a data frame of covariates covariates for the marginal regression models.")
  } else if(any(abd < 0) | any(abd >= 1)) {
    stop("ERROR: the entries of abd must be relative abundances between [0,1).")
  } else if(is.null(covars) == FALSE ) {
    if(all.equal(rownames(abd), rownames(covars)) != TRUE) {
      stop("ERROR: the rownames of abd and covars must match and be in the same order.")
    }
  }

  if(ncores > 1){
    message("Setting up paralleization")
    future::plan(future::multisession, workers = ncores)
  }

  message("Fitting marginal models.")
  if(is.null(covars)){
    mMLE = furrr::future_map(colnames(abd),
                             function(y, abd_mat) {
                               df = data.frame(y = abd_mat[,y])
                               out = marginalMLE(formula.mu = y ~ 1, formula.phi = ~ 1, formula.p = ~ 1, df = df)
                               return(out)
                             },
                             abd_mat = abd,
                             .progress = TRUE, .options = furrr::furrr_options(seed = T) )
  } else {
    mMLE = furrr::future_map(colnames(abd),
                             function(y, abd_mat, vars) {
                               df = cbind(y = abd_mat[,y], vars)
                               out = marginalMLE(formula.mu = y ~ ., formula.phi = ~ ., formula.p = ~ ., df = df)
                               return(out)
                             },
                             abd_mat = abd, vars = covars,
                             .progress = TRUE, .options = furrr::furrr_options(seed = T) )
  }

  names(mMLE) = colnames(abd)

  fitted.df = data.frame(p = I(lapply(mMLE, function(x) x$fitted$p)),
                         mu = I(lapply(mMLE, function(x) x$fitted$mu)),
                         phi = I(lapply(mMLE, function(x) x$fitted$phi)) )
  rownames(fitted.df) = colnames(abd)

  pairs = utils::combn(colnames(abd), 2)
  pairs = as.data.frame(t(pairs))
  colnames(pairs) = c("feature1", "feature2")

  pairs.data = furrr::future_map2(pairs$feature1, pairs$feature2,
                                  function(x, y, df) as.data.frame(df[,c(x, y)]),
                                  df = abd)

  pairs = cbind(pairs, data = I(pairs.data))

  pairs = merge(pairs, fitted.df, by.x = "feature1", by.y = "row.names", all.x = TRUE)
  colnames(pairs)[colnames(pairs) %in% c("p", "mu", "phi")] = c("p1", "mu1", "phi1")

  pairs = merge(pairs, fitted.df, by.x = "feature2", by.y = "row.names", all.x = TRUE)
  colnames(pairs)[colnames(pairs) %in% c("p", "mu", "phi")] = c("p2", "mu2", "phi2")

  pairs = pairs[,c("feature1", "feature2", "data", "p1", "mu1", "phi1", "p2", "mu2", "phi2")]

  pairs = pairs[with(pairs, order(feature1, feature2)),]


  message("Estimating copula dependence parameters.")
  theta = furrr::future_map_dbl(1:nrow(pairs),
                                function(dat, i){
                                  out = thetaTSMLE(x = dat$data[[i]],
                                                   lower = -200,
                                                   upper = 200,
                                                   p1 = dat$p1[[i]],
                                                   mu1 = dat$mu1[[i]],
                                                   phi1 = dat$phi1[[i]],
                                                   p2 = dat$p2[[i]],
                                                   mu2 = dat$mu2[[i]],
                                                   phi2 = dat$phi2[[i]]
                                  )

                                  return(out)
                                }, dat = pairs,
                                .progress = TRUE, .options = furrr::furrr_options(seed = T) )

  pairs = cbind(pairs, theta = theta)


  if(test == TRUE) {
    message("Performing likelihood ratio tests.")
    LR.ts = furrr::future_map_dbl(1:nrow(pairs),
                                  function(dat, i) {
                                    tryCatch({
                                      out = lrtTheta(response = dat$data[[i]],
                                                     p1 = dat$p1[[i]],
                                                     mu1 = dat$mu1[[i]],
                                                     phi1 = dat$phi1[[i]],
                                                     p2 = dat$p2[[i]],
                                                     mu2 = dat$mu2[[i]],
                                                     phi2 = dat$phi2[[i]],
                                                     theta = dat$theta[[i]],
                                                     thetaNull = 0)
                                    },
                                    error = function(e) {NA})
                                  }, dat = pairs,
                                  .progress = TRUE, .options = furrr::furrr_options(seed = T) )

    pairs = cbind(pairs, LR.ts = LR.ts)
  }

  return(pairs)
}
