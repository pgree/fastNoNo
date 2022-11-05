# fastNoNo

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/cmdstanr)](https://CRAN.R-project.org/package=fastNoNo)
[![R-CMD-check](https://github.com/pgree/fastNoNo/workflows/R-CMD-check/badge.svg)](https://github.com/pgree/fastNoNo/actions)
<!-- badges: end -->

An R package that fits fast approximations to Bayesian hierarchical linear
regression models using numerical linear algebra and low dimensional Gaussian
quadrature. This is an implementation of the method described in Greengard et
al. (2022).

* Philip Greengard, Jeremy Hoskins, Charles C. Margossian, Jonah Gabry, Andrew Gelman, and Aki Vehtari. (2022). Fast methods for posterior inference of two-group normal-normal models. [Bayesian Analysis](https://projecteuclid.org/journals/bayesian-analysis/advance-publication/Fast-Methods-for-Posterior-Inference-of-Two-Group-Normal-Normal/10.1214/22-BA1329.full)

### Installation

Installing from GitHub requires the necessary toolchain for compiling C++ code
via Rcpp (a pre-compiled version will be available from CRAN in the future). 

We recommend setting the compiler optimization option `-O3` in your `Makevars`
file in order to achieve the fastest speed possible when using the package. The
easiest way to do this is to use `usethis::edit_r_makevars()`, which will
automatically open the file you need to edit (or create it if it doesn't exist).
Alternatively, you can find or create the file yourself in your home directory
at `.R/Makevars`.

```r
# install.packages("usethis")
usethis::edit_r_makevars()

# this will open the Makevars file
# add the line CXXFLAGS += -O3
# then close the file
```

Then install the `fastNoNo` package from GitHub using
`remotes::install_github()`.

```r
# install.packages("remotes")
remotes::install_github("pgree/fastNoNo")
```


### Simple example
In the following code snippet, we fit a "mixed effects" model with the mtcars 
dataset that comes with R. We currently support two formats for specifying 
the mixed effects model and data. One is in the `lmer` style:
```
library(fastNoNo)

fit <- fit_mixed_formula(
  formula = mpg ~ wt + as.factor(gear) + (1|cyl),
  data = mtcars,
  sd_beta2 = 10,
  sd_sigma_y = 10,
  sd_sigma1 = 5,
  nnt = 40  # the number of quadrature nodes. increasing nnt will increase accuracy and runtime.
)
fit$beta1
fit$beta2
fit$sigma

# check accuracy of estimates
fit$errors

# refitting using more quadrature nodes improves accuracy
fit <- fit_mixed_formula(
  formula = mpg ~ wt + as.factor(gear) + (1 | cyl),
  data = mtcars,
  sd_beta2 = 10,
  sd_sigma_y = 10,
  sd_sigma1 = 5,
  nnt = 80
)
fit$errors
```

The other model/data-specification involves using data matrices
directly. In the following,
we use the stats::model.matrix function to construct the data matrices 
corresponding to the same data and model as the previous model. 
```
library(fastNoNo)

# suppose we want to fit the following model (in lme4 syntax)
# using the mtcars dataset that comes with R:
#   mpg ~ wt + as.factor(gear) + (1|cyl)
fit <- fit_mixed(
  y = mtcars$mpg,
  X1 = stats::model.matrix(~ 0 + as.factor(cyl), data = mtcars),
  X2 = stats::model.matrix(~ 1 + wt + as.factor(gear), data = mtcars),
  sd_sigma_y = 10,
  sd_sigma1 = 5,
  sd_beta2 = 10,
  nnt = 30
)
fit$beta1
fit$beta2
fit$sigma
```

### Example with simulated data
The following is a simulated data example using matrices and vectors 
as inputs for representing the data. 
```
library(fastNoNo)

n <- 10000
k1 <- 50
k2 <- 60

sigma_y <- abs(rnorm(1, 0, 1))
sigma1 <- abs(rnorm(1, 0, 1))
beta1 <- rnorm(k1, 0, sigma1)
beta2 <- rnorm(k2, 0, 1)

X1 <- matrix(rnorm(n * k1, 2, 3), ncol = k1)
X2 <- matrix(rnorm(n * k2, -1, 5), ncol = k2)
y <- rnorm(n, X1 %*% beta1 + X2 %*% beta2, sigma_y)

# Fit model
fit <- fit_mixed(
  y,
  X1,
  X2,
  sd_sigma_y = 1,
  sd_sigma1 = 1,
  sd_beta2 = rep(1, k2),
  nnt = 20
)
str(fit)

# Plot estimates of the betas vs "truth"
plot(fit$beta1$mean, beta1, pch = 20); abline(0, 1, col = "red")
plot(fit$beta2$mean, beta2, pch = 20); abline(0, 1, col = "red")
```
