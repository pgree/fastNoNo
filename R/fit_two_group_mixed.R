#' Fast two-group normal-normal model
#'
#' @description
#' Fits a fast and accurate approximation to the Bayesian two-group hierarchical
#' linear regression model
#'
#' \deqn{y ~ normal(X_1 \beta_1 + X_2 \beta_2, \sigma_y)}
#' \deqn{\beta_1 ~ normal(0, \sigma_1)}
#' \deqn{\beta_2 ~ normal(0, ss * I)}
#' \deqn{\sigma_y ~ normal+(0, sd_y)}
#' \deqn{\sigma_1 ~ normal+(0, sd_1)}
#'
#' where \eqn{ss} is a vector of positive numbers and \eqn{I} is the identity
#' matrix. The algorithm for computing the fit uses numerical linear algebra and
#' low dimensional Gaussian quadrature. See Greengard et al. (2021) for details.
#'
#' @export
#' @param y Outcome vector.
#' @param X1 Data matrix corresponding to group 1.
#' @param X2 Data matrix corresponding to group 2.
#' @param ss Vector of scale parameter priors corresponding to group 2.
#' @param sd_y Hyperprior on residual standard deviation.
#' @param sd1 Hyperprior on standard deviation of group 1.
#' @param nnt Number of quadrature nodes in \eqn{\theta}. See Greengard et al.
#'   (2021) for details.
#'
#' @return A named list with the following components:
#' * `beta1`: A data frame with two columns (`mean`, `sd`) containing the
#' posterior means and standard deviations for the vector \eqn{\beta_1} (the
#' coefficients on `X1`).
#' * `beta2`: A data frame with two columns (`mean`, `sd`)  containing the
#' posterior means and standard deviations for the vector \eqn{\beta_2} (the
#' coefficients on `X2`).
#' * `sigma`: A data frame with two columns (`mean`, `sd`) containing the
#' posterior means and standard deviations for \eqn{\sigma_y} and \eqn{\sigma_1}.
#' * `cov`: The posterior covariance matrix of coefficients \eqn{[\beta_1, \beta_2]}.
#' * `errors`: A data frame with two columns (`error_means`, `error_sds`)
#' containing the approximate accuracy of the posterior mean and standard
#' deviation estimates.
#' * `time`: The execution time of the algorithm in seconds. This only includes
#' the time taken to run the internal C++ code. To include the full time elapsed
#' when running the \R function use [system.time()].
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(1)
#' n <- 1000
#' k1 <- 50
#' k2 <- 60
#'
#' sigma_y <- 1
#' sigma1 <- 0.5
#' beta1 <- rnorm(k1, 0, sigma1)
#' beta2 <- rnorm(k2, 0, 1)
#'
#' X1 <- matrix(rnorm(n * k1, 2, 3), ncol = k1)
#' X2 <- matrix(rnorm(n * k2, -1, 5), ncol = k2)
#' y <- rnorm(n, X1 %*% beta1 + X2 %*% beta2, sigma_y)
#'
#' # Fit model
#' fit <- fit_two_group_mixed(y, X1, X2, ss = rep(1, k2), sd_y = 1, sd1 = 1, nnt = 20)
#' str(fit)
#'
#' # Plot estimates of the betas vs "truth"
#' plot(fit$beta1$mean, beta1); abline(0, 1, col = "red")
#' plot(fit$beta2$mean, beta2); abline(0, 1, col = "red")
#' }
#'
#' @references
#' Philip Greengard, Jeremy Hoskins, Charles C. Margossian, Jonah Gabry, Andrew
#' Gelman, and Aki Vehtari. (2022). Fast methods for posterior inference of
#' two-group normal-normal models. To appear,
#' [Bayesian Analysis](http://www.stat.columbia.edu/~gelman/research/published/two_group_fastnono.pdf)
#'
fit_two_group_mixed <- function(y, X1, X2, ss = rep(1, ncol(X2)), sd_y = 1, sd1 = 1, nnt = 10) {
  stopifnot(
    nrow(X1) == nrow(X2),
    length(y) == nrow(X1),
    length(ss) == ncol(X2),
    all(ss > 0),
    sd_y > 0,
    sd1 > 0,
    length(nnt) == 1,
    nnt >= 1
  )

  out1 <- run_two_group_mixed(y, X1, X2, ss, sd_y, sd1, nnt)
  out2 <- run_two_group_mixed(y, X1, X2, ss, sd_y, sd1, nnt = 2 * nnt)

  k1 <- ncol(X1)
  k2 <- ncol(X2)
  k <- k1 + k2

  # compute errors
  error_means <- out1$means - out2$means
  error_sds <- out1$sds - out2$sds
  errors <- data.frame(error_means, error_sds)
  rownames(errors) <- c(paste0("beta1_", 1:k1),
                        paste0("beta2_", 1:k2),
                        "sigma_y",
                        "sigma_beta1")

  # means for first group
  beta1 <- data.frame(out2$means[1:k1], out2$sds[1:k1])
  rownames(beta1) <- paste0("beta1_", 1:k1)
  colnames(beta1) <- c("mean", "sd")

  # means for second group
  beta2 <- data.frame(out2$means[(k1 + 1):k], out2$sds[(k1 + 1):k])
  rownames(beta2) <- paste0("beta2_", 1:k2)
  colnames(beta2) <- c("mean", "sd")

  # scale parameters
  sigma <- data.frame(out2$means[(k + 1):(k + 2)], out2$sds[(k + 1):(k + 2)])
  rownames(sigma) <- c("sigma_y", "sigma_beta1")
  colnames(sigma) <- c("mean", "sd")

  # posterior covariance of [beta1, beta2]
  cov <- matrix(data = out2$cov, nrow = k, ncol = k)

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

run_two_group_mixed <- function(y, X1, X2, ss, sd_y, sd1, nnt) {
  n <- length(y)
  k1 <- ncol(X1)
  k2 <- ncol(X2)
  fit <- mixed_2group_cpp(
     nnt = as.integer(nnt),  # number of quadrature in theta direction
     nn = as.integer(80),    # number of quadrature nodes in other directions
     n = as.integer(n),
     k1 = as.integer(k1),
     k2 = as.integer(k2),
     k = as.integer(k1+k2),
     a = cbind(X1, X2),
     y = y,
     ss = as.double(ss),
     sigy = as.double(sd_y),
     sig1 = as.double(sd1)
  )
  list(
    means = fit$means,
    sds = fit$sds,
    cov = fit$cov,
    time = fit$time
  )
}

