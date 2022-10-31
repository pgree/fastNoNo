#' Fast normal-normal mixed effects model
#'
#' @description
#' Fits a fast and accurate approximation to the Bayesian hierarchical
#' linear regression model
#'
#' \deqn{y ~ normal(X_1 \beta_1 + X_2 \beta_2, \sigma_y)}
#' \deqn{\beta_1 ~ normal(0, \sigma_1)}
#' \deqn{\beta_2 ~ normal(0, sd_beta2 * I)}
#' \deqn{\sigma_y ~ normal+(0, sd_sigma_y)}
#' \deqn{\sigma_1 ~ normal+(0, sd_1)}
#'
#' where \eqn{sd_beta2} is a vector of positive numbers and \eqn{I} is the identity
#' matrix. The algorithm for computing the fit uses numerical linear algebra and
#' low dimensional Gaussian quadrature. See Greengard et al. (2022) for details.
#'
#' @export
#' @param y (vector) The outcome variable.
#' @param X1 (matrix) The design matrix for the "random effects" part of the
#'   model. Must have `length(y)` rows.
#' @param X2 (matrix) The design matrix for the "fixed effects" part of the
#'   model. Must have `length(y)` rows.
#' @param sd_sigma_y (positive real) Scale parameter value for the prior on \eqn{\sigma_y}.
#' @param sd_sigma1 (positive real) Scale parameter value for the prior on \eqn{\sigma_1}.
#' @param sd_beta2 (positive reals) Scale parameter values for the prior on
#'   \eqn{\beta_2}. Must have either one element or `ncol(X2)` elements. In the
#'   former case the value is recycled.
#' @param nnt (positive integer) Number of quadrature nodes in \eqn{\theta}
#'   direction as described in Greengard et al. (2022).
#'
#' @return A named list with the following components:
#' * `beta1`: A data frame with two columns (`mean`, `sd`) containing the
#' posterior means and standard deviations for the vector \eqn{\beta_1}. If `X1`
#' has column names they are used as the row names for `beta1`, otherwise
#' generic names are used (`beta1_1`, `beta1_2`, etc.).
#' * `beta2`: A data frame with two columns (`mean`, `sd`)  containing the
#' posterior means and standard deviations for the vector \eqn{\beta_2}. If `X2`
#' has column names they are used as the row names for `beta2`, otherwise
#' generic names are used (`beta2_1`, `beta2_2`, etc.).
#' * `sigma`: A data frame with two columns (`mean`, `sd`) containing the
#' posterior means and standard deviations for \eqn{\sigma_y} and \eqn{\sigma_1}.
#' * `cov`: The posterior covariance matrix of coefficients \eqn{[\beta_1, \beta_2]}.
#' * `errors`: A data frame with two columns (`error_mean`, `error_sd`)
#' containing the approximate accuracy of the posterior mean and standard
#' deviation estimates.
#' * `time`: The execution time of the algorithm in seconds. This only includes
#' the time taken to run the internal C++ code. To include the full time elapsed
#' when running the \R function see [system.time()].
#'
#' @examples
#' ### Example with simulated data
#' set.seed(1)
#' n <- 1000
#' k1 <- 50
#' k2 <- 60
#'
#' sigma_y <- abs(rnorm(1, 0, 1))
#' sigma1 <- abs(rnorm(1, 0, 1))
#' beta1 <- rnorm(k1, 0, sigma1)
#' beta2 <- rnorm(k2, 0, 1)
#'
#' X1 <- matrix(rnorm(n * k1, 2, 3), ncol = k1)
#' X2 <- matrix(rnorm(n * k2, -1, 5), ncol = k2)
#' y <- rnorm(n, X1 %*% beta1 + X2 %*% beta2, sigma_y)
#'
#' # Fit model
#' fit <- fit_two_group_mixed(y, X1, X2, sd_beta2 = rep(1, k2), sd_sigma_y = 1, sd_sigma1 = 1, nnt = 20)
#' str(fit)
#'
#' # Plot estimates of the betas vs "truth"
#' plot(fit$beta1$mean, beta1, pch = 20); abline(0, 1, col = "red")
#' plot(fit$beta2$mean, beta2, pch = 20); abline(0, 1, col = "red")
#'
#'
#' ### Example of varying intercept model using mtcars dataset
#'
#' # suppose we want to fit the following model (in lme4 syntax)
#' # using the mtcars dataset that comes with R:
#' #   mpg ~ wt + as.factor(gear) + (1|cyl)
#'
#' fit <- fit_two_group_mixed(
#'   y = mtcars$mpg,
#'   X1 = stats::model.matrix(~ 0 + as.factor(cyl), data = mtcars),
#'   X2 = stats::model.matrix(~ wt + as.factor(gear), data = mtcars),
#'   sd_sigma_y = 10,
#'   sd_sigma1 = 5,
#'   sd_beta2 = 10,
#'   nnt = 30
#' )
#' fit$beta1
#' fit$beta2
#' fit$sigma
#'
#' @references
#' Philip Greengard, Jeremy Hoskins, Charles C. Margossian, Jonah Gabry, Andrew
#' Gelman, and Aki Vehtari. (2022). Fast methods for posterior inference of
#' two-group normal-normal models. To appear,
#' [Bayesian Analysis](http://www.stat.columbia.edu/~gelman/research/published/two_group_fastnono.pdf)
#'
fit_two_group_mixed <- function(y, X1, X2, sd_sigma_y = 1, sd_sigma1 = 1, sd_beta2 = rep(1, ncol(X2)), nnt = 10) {
  stopifnot(
    !anyNA(y),
    !anyNA(X1),
    !anyNA(X2),
    nrow(X1) == nrow(X2),
    length(y) == nrow(X1),
    length(sd_beta2) == 1 || length(sd_beta2) == ncol(X2),
    length(sd_sigma_y) == 1,
    length(sd_sigma1) == 1,
    length(nnt) == 1,
    all(sd_beta2 > 0),
    sd_sigma_y > 0,
    sd_sigma1 > 0,
    nnt > 0,
    nnt == as.integer(nnt)
  )
  if (length(sd_beta2) == 1) {
    sd_beta2 <- rep(sd_beta2, ncol(X2))
  }

  out1 <- run_two_group_mixed(y, X1, X2, sd_beta2, sd_sigma_y, sd_sigma1, nnt)
  out2 <- run_two_group_mixed(y, X1, X2, sd_beta2, sd_sigma_y, sd_sigma1, nnt = 2 * nnt)

  k1 <- ncol(X1)
  k2 <- ncol(X2)
  k <- k1 + k2

  beta1 <- data.frame(out2$means[1:k1], out2$sds[1:k1])
  rownames(beta1) <- if (!is.null(colnames(X1))) colnames(X1) else paste0("beta1_", 1:k1)
  colnames(beta1) <- c("mean", "sd")

  beta2 <- data.frame(out2$means[(k1 + 1):k], out2$sds[(k1 + 1):k])
  rownames(beta2) <- if (!is.null(colnames(X2))) colnames(X2) else paste0("beta2_", 1:k2)
  colnames(beta2) <- c("mean", "sd")

  sigma <- data.frame(out2$means[(k + 1):(k + 2)], out2$sds[(k + 1):(k + 2)])
  rownames(sigma) <- c("sigma_y", "sigma_beta1")
  colnames(sigma) <- c("mean", "sd")

  errors <- data.frame(out1$means - out2$means, out1$sds - out2$sds)
  rownames(errors) <- c(rownames(beta1), rownames(beta2), "sigma_y", "sigma_beta1")
  colnames(errors) <- c("error_mean", "error_sd")

  cov <- matrix(data = out2$cov, nrow = k, ncol = k)
  rownames(cov) <- c(rownames(beta1), rownames(beta2))
  colnames(cov) <- rownames(cov)

  list(
    beta1 = beta1,
    beta2 = beta2,
    sigma = sigma,
    cov = cov,
    errors = errors,
    time = out1$time
  )
}


# internal ----------------------------------------------------------------

# run the c++ code for the fastNoNo algorithm
run_two_group_mixed <- function(y, X1, X2, sd_beta2, sd_sigma_y, sd_sigma1, nnt) {
  mixed_2group_cpp(
     nnt = as.integer(nnt),  # number of quadrature in theta direction
     nn = as.integer(80),    # number of quadrature nodes in other directions
     n = length(y),
     k1 = ncol(X1),
     k2 = ncol(X2),
     k = ncol(X1) + ncol(X2),
     a = cbind(X1, X2),
     y = y,
     sd_beta2 = as.double(sd_beta2),
     sigy = as.double(sd_sigma_y),
     sig1 = as.double(sd_sigma1)
  )
}

