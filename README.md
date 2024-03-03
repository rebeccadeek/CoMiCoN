
# `CoMiCoN`: Copulas with Mixture Margins for Covariation Networks

<!-- badges: start -->
<!-- badges: end -->

CoMiCoN is an R package for microbial covariation and network analysis. It implements the copula model described in [Inference of microbial covariation networks using copula models with mixture margins](https://academic.oup.com/bioinformatics/article/39/7/btad413/7209520) by Deek and Li (2023), now in published in *Bioinformatics*.

## About

CoMiCoN uses a mixture margin bivariate copula model for a pair of zero-inflated microbial relative abundances ($X_i$, $X_j$): 
$F(x_i, x_j; \boldsymbol{\gamma}_i, \boldsymbol{\gamma}_j, \theta_{ij}) = C(F_i(x_i; \boldsymbol{\gamma}_i), F_j(x_j; \boldsymbol{\gamma}_j); \theta_{ij})$

The copula dependence parameter $\theta$ captures the covariation for the pair of microbes and can be used to build covariation networks. The package provides implementation of a two-stage maximum likelihood estimation (tsMLE) procedure and a two-stage likelihood ratio test (tsLRT) for the CoMiCoN model.

## Installation

You can install the latest version of `CoMiCoN` from GitHub with:

``` r
install.packages("devtools")
devtools::install_github("rebeccadeek/CoMiCoN")
```

## Documentation and Examples

Help documentation for the `CoMiCoN` package is available in R. After installing the package from GitHub via `devtools` and loading it with `library()` use `?` to access the documentation of functions in the package. E.g.

``` r
?dzib
```

The main function, `comicon` takes two data frames:

1. `abd`: a $n \times p$ data frame of relative abundance, with samples as rows and microbes as columns.
2. `covars`: a $n \times q$ data frame of covariates, with samples as rows in the same order as abd.

The function has two other inputs, `test`: a logical indicating if the tsLRT should be run and `ncores`: the number of cores for parallelization. The default `ncores` is one core, implying no parallel computing.

## Contact

To report any bugs, issues, or suggestions please contact the maintainer Rebecca Deek via [email](mailto:rdeek@pitt.edu).
