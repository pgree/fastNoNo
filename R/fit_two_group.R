#' Fast two-group normal-normal model
#'
#' @description
#' This function fits a fast approximation to the Bayesian two-group
#' hierarchical linear regression model
#'
#' \deqn{y ~ normal(X_1 \beta_1 + X_2 \beta_2, \sigma_y)}
#' \deqn{\beta_1 ~ normal(0, \sigma_1)}
#' \deqn{\beta_2 ~ normal(0, \sigma_2)}
#' \deqn{\sigma_y ~ normal+(0, sd_y)}
#' \deqn{\sigma_1 ~ normal+(0, sd_1)}
#' \deqn{\sigma_2 ~ normal+(0, sd_2)}
#'
#' The algorithm for computing the fit uses numerical linear algebra and
#' low dimensional Gaussian quadrature. See Greengard et al. (2021) for details.
#'
#' @export
#' @param y Outcome vector.
#' @param X1 Data matrix corresponding to group 1.
#' @param X2 Data matrix corresponding to group 2.
#' @param sd_y Prior on residual standard deviation.
#' @param sd_1 Hyperprior on standard deviation of group 1.
#' @param sd_2 Hyperprior on standard deviation of group 2.
#' @param nnt Number of quadrature nodes in \eqn{\theta}. See Greengard et al.
#'   (2021) for details.
#'
#' @return A named list with the following components:
#' * `beta_1`: A data frame of posterior means and standard deviations for
#' the vector \eqn{\beta_1}, the coefficients on `X1`.
#' * `beta_2`: A data frame of posterior means and standard deviations for
#' the vector \eqn{\beta_2}, the coefficients on `X2`.
#' * `sigma`: A data frame with posterior means and standard deviations for
#' \eqn{\sigma_y}, \eqn{\sigma_1}, and \eqn{\sigma_2}.
#' * `errors`: A data frame with approximate accuracy of the posterior
#' mean and standard deviation estimates.
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' n <- 10000
#' k_1 <- 50
#' k_2 <- 60
#'
#' sigma_y <- 1
#' sigma_1 <- 0.5
#' sigma_2 <- 2
#' beta_1 <- rnorm(k_1, 0, sigma_1)
#' beta_2 <- rnorm(k_2, 0, sigma_2)
#'
#' X_1 <- matrix(rnorm(n * k_1, 2, 3), ncol = k_1)
#' X_2 <- matrix(rnorm(n * k_2, -1, 5), ncol = k_2)
#' y <- rnorm(n, X_1 %*% beta_1 + X_2 %*% beta_2, sigma_y)
#'
#' # hyperpriors on scale parameters
#' sd_1 <- 1
#' sd_2 <- 2
#' sd_y <- 2
#'
#' # Fit model
#' fit <- fit_two_group(y, X_1, X_2, sd_y, sd_1, sd_2)
#' str(fit)
#'
#' # Plot estimates of the betas vs "truth"
#' plot(fit$beta_1$mean, beta_1); abline(0, 1, col = "red")
#' plot(fit$beta_2$mean, beta_2); abline(0, 1, col = "red")
#' }
#'
#' @references
#' Philip Greengard, Jeremy Hoskins, Charles C. Margossian, Andrew Gelman, and Aki Vehtari.
#' (2021). Fast methods for posterior inference of two-group normal-normal models.
#' [preprint arXiv:2110.03055](https://arxiv.org/abs/2110.03055)
#'
#' @useDynLib fastNoNo dense_eval
fit_two_group <- function(y, X1, X2, sd_y, sd_1 = 1, sd_2 = 1, nnt = 10) {
  stopifnot(nrow(X1) == nrow(X2), length(y) == nrow(X1), sd_y > 0,
            sd_1 > 0, sd_2 > 0, length(nnt) == 1, nnt >= 1)

  out1 <- run_two_group(y, X1, X2, sd_y, sd_1, sd_2, nnt)
  out2 <- run_two_group(y, X1, X2, sd_y, sd_1, sd_2, nnt = 2 * nnt)

  k1 <- ncol(X1)
  k2 <- ncol(X2)
  k <- k1 + k2

  # compute errors
  error_means <- out1$means - out2$means
  error_sds <- out1$sds - out2$sds
  errors <- data.frame(error_means, error_sds)
  rownames(errors) <- c(paste0("beta_1_", 1:k1), paste0("beta_2_", 1:k2),
                        "sigma_y", "sigma_1", "sigma_2")

  # means for first group
  beta_1 <- data.frame(out2$means[1:k1], out2$sds[1:k1])
  rownames(beta_1) <- paste0("beta_1_", 1:k1)
  colnames(beta_1) <- c("mean", "sd")

  # means for second group
  beta_2 <- data.frame(out2$means[(k1 + 1):k], out2$sds[(k1 + 1):k])
  rownames(beta_2) <- paste0("beta_2_", 1:k2)
  colnames(beta_2) <- c("mean", "sd")

  # scale parameters
  sigma <- data.frame(out2$means[(k + 1):(k + 3)], out2$sds[(k + 1):(k + 3)])
  rownames(sigma) <- c("sigma_y", "sigma_1", "sigma_2")
  colnames(sigma) <- c("mean", "sd")

  list(
    beta_1 = beta_1,
    beta_2 = beta_2,
    sigma = sigma,
    errors = errors
  )
}


# internal ----------------------------------------------------------------

run_two_group <- function(y, X1, X2, sd_y, sd_1, sd_2, nnt) {
  # extract parameters from inputs
  n <- length(y)
  k1 <- ncol(X1)
  k2 <- ncol(X2)

  fit <- .Fortran(
    "dense_eval",
    nnt = as.integer(nnt),  # number of quadrature in theta direction
    nn = as.integer(80),    # number of quadrature nodes in other directions
    n = as.integer(n),
    k1 = as.integer(k1),
    k2 = as.integer(k2),
    X = cbind(X1, X2),
    y = y,
    sigs = as.double(c(sd_y, sd_1, sd_2)),
    # these are dummy objects for fortran to use for the results
    means = as.double(rep(-99, k1 + k2 + 3)),
    sds = as.double(rep(-99, k1 + k2 + 3))
  )

  list(means = fit$means, sds = fit$sds)
}
