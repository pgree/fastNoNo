# fastNoNo

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/cmdstanr)](https://CRAN.R-project.org/package=fastNoNo)
[![R-CMD-check](https://github.com/pgree/fastNoNo/workflows/R-CMD-check/badge.svg)](https://github.com/pgree/fastNoNo/actions)
<!-- badges: end -->

An R package that fits fast approximations to Bayesian hierarchical linear
regression models using numerical linear algebra and low dimensional Gaussian
quadrature. For more details see the following paper:

* Philip Greengard, Jeremy Hoskins, Charles C. Margossian, Jonah Gabry, Andrew Gelman, and Aki Vehtari. (2022). Fast methods for posterior inference of two-group normal-normal models. To appear, [Bayesian Analysis](http://www.stat.columbia.edu/~gelman/research/published/two_group_fastnono.pdf)

### Installation

Installing from GitHub requires the necessary toolchain for compiling C++ code
via Rcpp (a pre-compiled version will be available from CRAN in the future). 

We recommend setting the compiler optimization option `-O3` in your `Makevars`
file in order to achieve the fastest speed possible when using the package. The
easiest way to do this is to use the function `edit_r_makevars()` from the
`usethis` package, which will automatically open the file you need to edit (or
create it if it doesn't exist). Alternatively you can find or create the file
yourself in your home directory at `.R/Makevars`.

```r
# install.packages("usethis")
usethis::edit_r_makevars()

# this will open the Makevars file
# add the line CXXFLAGS += -O3
# then close the file
```

Then install the `fastNoNo` package from GitHub using the `remotes` package. 

```r
# install.packages("remotes")
remotes::install_github("pgree/fastNoNo", ref = "cpp")
```
