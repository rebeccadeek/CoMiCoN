# load packages
library(tidyverse)
library(furrr)
library(gamlss)

# set up parallelization
nCores = future::availableCores() - 2
future::plan(future::multisession, workers = nCores)

# path to data
# path = "/path/to/RDS-file/"

# read data
g.relAbun20.meta.h = readRDS(paste0(path, "agp-cleaned-relAbun-metadata.RDS"))

# Estimation of marginal MLEs
# Function to calculate marginal zero-inflated beta MLEs from a ZIB regression 
# model returns a list with 3 components:
# 1. Regression coefficients; 2. Fitted p/mu/phi; 3. Fitted shape parameters alpha/beta

source("01b_20210304_twostageMLEFunction.R")

mmmleTblFunc = function(y, x.mu, x.phi, x.p, df){
  formula.mu = as.formula(paste0(y, paste0(x.mu, collapse = " "))) # paste "y" "~ x" together to make formula
  x.phi = as.formula(x.phi); x.p = as.formula(x.p)                 # convert to type=formula
  
  out = marginalMLE(formula.mu, x.phi, x.p, df)                    # marginal MLEs
  
  return(out)
}

# Estimate MLEs with covariate adjustment for: 
# 1. age; 2. BMI; 3. Antibiotic use

g.relAbun20.meta.h.mmle.C = tibble::tibble(genus = colnames(g.relAbun20.meta.h)[-c(1:4)]) %>% 
  dplyr::mutate(., mmle =                                                                                        # map over y = genus
                  furrr::future_map(genus, mmmleTblFunc, x.mu = ~ age + bmi + antibiotic_select,                 # covars for p 
                                    x.phi = ~ 1, x.p = ~ age + bmi + antibiotic_select, df = g.relAbun20.meta.h, # and mu
                                    .progress = TRUE, .options = furrr::furrr_options(seed = T)) )

# Taxon-taxon pairs
# Find all possible pairs of taxa

pw.genus = combn(colnames(g.relAbun20.meta.h)[-c(1:4)], 2)
rownames(pw.genus) = c("g1", "g2")
pw.genus = tibble::as_tibble(t(pw.genus))

# Estimation of dependence parameters
# Tidyverse function to merge marginal MLEs with taxon pairs + estimate dependence parameter (theta)
pwAssocFunc = function(pw, dat, margML) {
  # LRT try function
  lrtTheta.try = function(data, eta.shape) {
    source("02a_20210304_lrtFunction.R")
    
    tryCatch({
      lrtTheta(data, eta.shape)
    }, error = function(e) {
      NA
    } )
  }
  
  pw.assoc = pw %>% 
    dplyr::mutate(., data = furrr::future_map2(g1, g2, function(x, y, df) as.data.frame(df[,c(x, y)]), df = dat)) %>%  # pw relAb
    # microbe 1
    dplyr::left_join(., margML, by = c("g1" = "genus")) %>%          # merge pw relAb w margML1
    dplyr::mutate(., mmle = purrr::map(mmle, function(x) {
      names(x[["coefficients"]]) = paste0(names(x[["coefficients"]]), 1) # rename regression
      names(x[["fitted"]]) = paste0(names(x[["fitted"]]), 1)         # and fitted
      names(x[["shape"]]) = paste0(names(x[["shape"]]), 1)           # and shape components
      names(x) = paste0(names(x), 1)                                 # for 1st microbe
      return(x) }) ) %>%
    tidyr::unnest_wider(mmle) %>%                                    # unnest mMLE
    # microbe 2
    dplyr::left_join(., margML, by = c("g2" = "genus")) %>%          # merge pw relAb w margML2
    dplyr::mutate(., mmle = purrr::map(mmle, function(x) {
      names(x[["coefficients"]]) = paste0(names(x[["coefficients"]]), 2) # rename regression
      names(x[["fitted"]]) = paste0(names(x[["fitted"]]), 2)         # and fitted
      names(x[["shape"]]) = paste0(names(x[["shape"]]), 2)           # and shape components
      names(x) = paste0(names(x), 2)                                 # for 2nd microbe
      return(x)}) ) %>% 
    tidyr::unnest_wider(mmle) %>%                                    # unnest mMLE
    tidyr::unnest_wider(shape1) %>% tidyr::unnest_wider(shape2) %>%  # unnest shape
    dplyr::mutate(., theta = furrr::future_pmap_dbl(                                       # theta TSML
      dplyr::select(., c(data,p1,alpha1,beta1,p2,alpha2,beta2)) %>% dplyr::rename(x=data), # select and rename
      optTheta, lower = -200, upper = 200, .progress = TRUE,                               # theta opt
      .options = furrr::furrr_options(seed = T)) ) %>%                                    # pb & seed
    dplyr::mutate(., eta.shape = purrr::pmap(select(., p1,alpha1,beta1,p2,alpha2,beta2,theta), data.frame) ) %>%    # eta df
    dplyr::mutate(., LR.ts = furrr::future_map2_dbl(data, eta.shape, lrtTheta.try,                                  # LR ts
                                                    .progress = TRUE, .options = furrr::furrr_options(seed = T)) ) # pb & seed
  
  return(pw.assoc)
}

# Estimate dependence parameters using marginal MLEs calculated using above function
pw.genus.assoc.C = pwAssocFunc(pw = pw.genus, dat = g.relAbun20.meta.h, margML = g.relAbun20.meta.h.mmle.C)