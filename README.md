# fastNoNo

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/cmdstanr)](https://CRAN.R-project.org/package=fastNoNo)
[![R-CMD-check](https://github.com/pgree/fastNoNo/workflows/R-CMD-check/badge.svg)](https://github.com/pgree/fastNoNo/actions)
<!-- badges: end -->

**Note: this package is a work in progress and is not yet officially supported**

An R package that fits fast approximations to Bayesian hierarchical linear regression models 
using numerical linear algebra and low dimensional Gaussian quadrature.

### Installation

Installing from GitHub requires the necessary toolchain for compiling C++ code
via Rcpp and RcppEigen. A pre-compiled version will be available from CRAN in
the future.

```r
# install.packages("remotes")
remotes::install_github("pgree/fastNoNo")
```


