#' Fast two-group normal-normal model
#'
#' @description
#' This function fits a fast approximation to the
#' Bayesian two-group hierarchical linear regression model
#'
#' \deqn{y ~ normal(X_1 \beta_1 + X_2 \beta_2, \sigma_y)}
#' \deqn{\beta_1 ~ normal(0, \sigma_1)}
#' \deqn{\beta_2 ~ normal(0, \sigma_{beta_2})}
#' \deqn{\sigma_beta_1 ~ normal+(0, 1)}
#' \deqn{\sigma_y ~ normal+(0, 1)}
#'
#' using numerical linear algebra and low dimensional Gaussian quadrature.
#'
#' @export
#' @param y Outcome vector.
#' @param X1 Data matrix corresponding to group 1.
#' @param X2 Data matrix corresponding to group 2.
#' @param sigma_beta_2 Data vector of scale parameter priors corresponding to group 2.
#' @param nnt Number of quadrature nodes in \eqn{\theta}. See Greengard et al.
#'   for details.
#'
#' @return A named list with the following components:
#' * `beta_1`: A data frame of posterior means and standard deviations for
#' the vector \eqn{\beta_1}, the coefficients on `X1`.
#' * `beta_2`: A data frame of posterior means and standard deviations for
#' the vector \eqn{\beta_2}, the coefficients on `X2`.
#' * `sigma`: A data frame with posterior means and standard deviations for
#' \eqn{\sigma_y} and \eqn{\sigma_beta_1}.
#' * `cov`: The posterior covariance matrix of `beta_1` and `beta_2`.
#' * `errors_means`: A data frame with approximate accuracy of the posterior
#' mean and standard deviation estimates
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(1)
#' n <- 1000
#' k_1 <- 50
#' k_2 <- 60
#'
#' sigma_y <- 1
#' sigma_beta_1 <- 0.5
#' beta_1 <- rnorm(k_1, 0, sigma_beta_1)
#' beta_2 <- rnorm(k_2, 0, 1)
#' sigma_beta_2 <- rep(1.0, k_2)
#'
#' X_1 <- matrix(rnorm(n * k_1, 2, 3), ncol = k_1)
#' X_2 <- matrix(rnorm(n * k_2, -1, 5), ncol = k_2)
#' y <- rnorm(n, X_1 %*% beta_1 + X_2 %*% beta_2, sigma_y)
#'
#' # Fit model
#' fit <- fit_two_group_mixed(y, X_1, X_2, sigma_beta_2, nnt=20)
#' str(fit)
#'
#' # Plot estimates vs truth
#' plot(fit$beta_1$mean, beta_1)
#' plot(fit$beta_2$mean, beta_2)
#' plot(fit$sigma$mean, c(sigma_y, sigma_beta_1))
#' }
#'
#' @useDynLib fastNoNo dense_eval
fit_two_group_mixed <- function(y, X1, X2, ss, nnt = 10) {
  stopifnot(nrow(X1) == nrow(X2), length(y) == nrow(X1),
            length(nnt) == 1)

  out1 <- run_two_group_mixed(y, X1, X2, ss, nnt)
  out2 <- run_two_group_mixed(y, X1, X2, ss, nnt = 2*nnt)

  k1 <- ncol(X1)
  k2 <- ncol(X2)
  k <- k1+k2

  # compute errors
  error_means <- out1$means - out2$means
  error_sds <- out1$sds - out2$sds
  errors <- data.frame(error_means, error_sds)
  rownames(errors) <- c(paste0("beta_1_", 1:k1), paste0("beta_2_", 1:k2),
                        "sigma_y", "sigma_beta_1")

  # means for first group
  beta_1 <- data.frame(out2$means[1:k1], out2$sds[1:k1])
  rownames(beta_1) <- paste0("beta_1_", 1:k1)
  colnames(beta_1) <- c("mean", "sd")

  # means for second group
  beta_2 <- data.frame(out2$means[(k1+1):k], out2$sds[(k1+1):k])
  rownames(beta_2) <- paste0("beta_2_", 1:k2)
  colnames(beta_2) <- c("mean", "sd")

  # scale parameters
  sigma <- data.frame(out2$means[(k+1):(k+2)], out2$sds[(k+1):(k+2)])
  rownames(sigma) <- c("sigma_y", "sigma_beta_1")
  colnames(sigma) <- c("mean", "sd")

  # posterior covariance of beta_1 and beta_2
  cov <- matrix(data=out2$cov, nrow=k, ncol=k)

  list(
    beta_1 = beta_1,
    beta_2 = beta_2,
    sigma = sigma,
    cov = cov,
    errors = errors
  )
}


# internal ----------------------------------------------------------------

run_two_group_mixed <- function(y, X1, X2, ss, nnt) {
  # extract parameters from inputs
  n <- length(y)
  k1 <- ncol(X1)
  k2 <- ncol(X2)

  fit <- .Fortran(
    "mixed_2group",
    nnt = as.integer(nnt),  # number of quadrature in theta direction
    nn = as.integer(80),    # number of quadrature nodes in other directions
    n = n,
    k1 = as.integer(k1),
    k2 = as.integer(k2),
    k = as.integer(k1+k2),
    X = cbind(X1, X2),
    y = y,
    # these are dummy objects for fortran to use for the results
    ss = as.double(ss),
    means = as.double(rep(-99, k1+k2+2)),
    dsum = 0.0,
    sds = as.double(rep(-99, k1+k2+2)),
    cov = as.double(rep(-99, (k1+k2)**2))
  )

  list(means=fit$means, sds=fit$sds, cov=fit$cov)
}
