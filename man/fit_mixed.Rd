% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_mixed.R
\name{fit_mixed}
\alias{fit_mixed}
\title{Fast normal-normal mixed effects model}
\usage{
fit_mixed(y, X1, X2, ..., sd_sigma_y, sd_sigma1, sd_beta2, nnt = 20)
}
\arguments{
\item{y}{(vector) The outcome variable.}

\item{X1}{(matrix) The design matrix for the "random effects" part of the
model. Must have \code{length(y)} rows.}

\item{X2}{(matrix) The design matrix for the "fixed effects" part of the
model. Must have \code{length(y)} rows.}

\item{...}{Currently for internal use only.}

\item{sd_sigma_y}{(positive real) Scale parameter value for the prior on
\eqn{\sigma_y}.}

\item{sd_sigma1}{(positive real) Scale parameter value for the prior on
\eqn{\sigma1}.}

\item{sd_beta2}{(positive reals) Scale parameter values for the priors on
\eqn{\beta2} ("fixed effects" coefficients). Must have either one element
or \code{ncol(X2)} elements. In the former case the value is recycled.}

\item{nnt}{(positive integer) Number of quadrature nodes in \eqn{\theta}
direction as described in Greengard et al. (2022). We recommend increasing
\code{nnt} to improve the accuracy of the estimates if the errors are too large
(see \code{errors} slot in returned fitted model object).}
}
\value{
A named list with the following components:
\itemize{
\item \code{beta1}: A data frame with two columns (\code{mean}, \code{sd}) containing the
posterior means and standard deviations for the vector \eqn{\beta1} (the
"random effects").
\item \code{beta2}: A data frame with two columns (\code{mean}, \code{sd})  containing the
posterior means and standard deviations for the vector \eqn{\beta2} (the
"fixed effects").
\item \code{sigma}: A data frame with two columns (\code{mean}, \code{sd}) containing the
posterior means and standard deviations for \eqn{\sigma_y} (the residual
standard deviation) and \eqn{\sigma1} (the standard deviation of
\eqn{\beta1}).
\item \code{cov}: The posterior covariance matrix of \eqn{[\beta1, \beta2]}.
\item \code{errors}: A data frame with two columns (\code{error_mean}, \code{error_sd})
containing the approximate accuracy of the posterior mean and standard
deviation estimates.
\item \code{time}: The execution time of the algorithm in seconds. This only includes
the time taken to run the internal C++ code. To include the full time elapsed
when running the \R function see \code{\link[=system.time]{system.time()}}.
}
}
\description{
Fits a fast and accurate approximation to the Bayesian linear regression
model

\deqn{y ~ normal(X_1 \beta_1 + X_2 \beta_2, \sigma_y)}
\deqn{\beta1 ~ normal(0, \sigma1)}
\deqn{\beta2 ~ normal(0, sd_\beta2 * I)}
\deqn{\sigma_y ~ normal+(0, sd_\sigma_y)}
\deqn{\sigma1 ~ normal+(0, sd_\sigma1)}

where \eqn{sd_\beta2} is a vector of positive numbers and \eqn{I} is the
identity matrix. The algorithm for computing the fit uses numerical linear
algebra and low dimensional Gaussian quadrature. See Greengard et al. (2022)
for details.
}
\examples{
### Example with simulated data
set.seed(1)
n <- 10000
k1 <- 50
k2 <- 60

sigma_y <- abs(rnorm(1, 0, 1))
sigma1 <- abs(rnorm(1, 0, 1))
beta1 <- rnorm(k1, 0, sigma1)
beta2 <- rnorm(k2, 0, 1)

X1 <- matrix(rnorm(n * k1, 2, 3), ncol = k1)
X2 <- matrix(rnorm(n * k2, -1, 5), ncol = k2)
y <- rnorm(n, X1 \%*\% beta1 + X2 \%*\% beta2, sigma_y)

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


### Example of varying intercept model using mtcars dataset

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

}
\references{
Philip Greengard, Jeremy Hoskins, Charles C. Margossian, Jonah Gabry, Andrew
Gelman, and Aki Vehtari. (2022). Fast methods for posterior inference of
two-group normal-normal models. To appear,
\href{http://www.stat.columbia.edu/~gelman/research/published/two_group_fastnono.pdf}{Bayesian Analysis}
}
\seealso{
\code{\link[=fit_mixed_formula]{fit_mixed_formula()}}
}
