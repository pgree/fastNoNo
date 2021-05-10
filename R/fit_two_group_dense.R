#' Fast two-group normal-normal model
#'
#' @description
#' This function fits a fast approximation to the
#' Bayesian two-group hierarchical linear regression model
#'
#' \deqn{y ~ normal(X_1 \beta_1 + X_2 \beta_2, \sigma_y)}
#' \deqn{\beta_1 ~ normal(0, \sigma_1)}
#' \deqn{\beta_2 ~ normal(0, \sigma_2)}
#' \deqn{\sigma_1 ~ normal+(0, 1)}
#' \deqn{\sigma_2 ~ normal+(0, 1)}
#' \deqn{\sigma_y ~ normal+(0, 1)}
#'
#' using numerical linear algebra and low dimensional Gaussian quadrature.
#'
#' **NOTE:** Currently the algorithm is implemented in Fortran and must be
#' compiled the first time the function is called in an R session. Subsequent
#' calls to the function in the same R session do not require recompilation.
#' In future versions of the package the Fortran code will come pre-compiled.
#'
#' @export
#' @param y Outcome vector.
#' @param X1 Data matrix corresponding to group 1.
#' @param X2 Data matrix corresponding to group 2.
#'
#' @return A named list with the following components:
#' * `beta_1`: A data frame of posterior means and standard deviations for
#' the vector \eqn{\beta_1}, the coefficients on `X1`.
#' * `beta_2`: A data frame of posterior means and standard deviations for
#' the vector \eqn{\beta_2}, the coefficients on `X2`.
#' * `sigma`: A data frame with posterior means and standard deviations for
#' \eqn{\sigma_y}, \eqn{\sigma_1}, and \eqn{\sigma_2}.
#' * `errors_means`: A data frame with approximate accuracy of the posterior
#' mean and standard deviation estimates
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
#' # Fit model
#' fit <- fit_two_group_dense(y, X_1, X_2)
#' str(fit)
#'
#' # Plot estimates vs truth
#' plot(fit$beta_1$mean, beta_1)
#' plot(fit$beta_2$mean, beta_2)
#' plot(fit$sigma$mean, c(sigma_y, sigma_1, sigma_2))
#' }
#'
#' @useDynLib fastNoNo dense_eval
fit_two_group_dense <- function(y, X1, X2) {
  stopifnot(nrow(X1) == nrow(X2), length(y) == nrow(X1))

  # nodes in theta direction
  nnt <- 10
  out1 <- run_two_group_dense(y, X1, X2, nnt)

  nnt2 <- 2*nnt
  out2 <- run_two_group_dense(y, X1, X2, nnt2)

  # process output
  k1 <- ncol(X1)
  k2 <- ncol(X2)
  k <- k1+k2

  # compute errors
  error_means <- out1$means - out2$means
  error_stds <- out1$stds - out2$stds
  errors <- data.frame(error_means, error_stds)
  print(errors)
  print(c(paste0("beta_1_", 1:k1), paste0("beta_2_", 1:k2)))
  rownames(errors) <- c(paste0("beta_1_", 1:k1), paste0("beta_2_", 1:k2),
                        "sigma_y", "sigma_1", "sigma_2")

  # means for first group
  beta_1 <- data.frame(out2$means[1:k1], out2$stds[1:k1])
  rownames(beta_1) <- paste0("beta_1_", 1:k1)
  colnames(beta_1) <- c("mean", "std")

  # means for second group
  beta_2 <- data.frame(out2$means[(k1+1):k], out2$stds[(k1+1):k])
  rownames(beta_2) <- paste0("beta_2_", 1:k2)
  colnames(beta_2) <- c("mean", "std")

  # scale parameters
  sigma <- data.frame(out2$means[(k+1):(k+3)], out2$stds[(k+1):(k+3)])
  rownames(sigma) <- c("sigma_y", "sigma_1", "sigma_2")
  colnames(sigma) <- c("mean", "std")

  list(
    beta_1 = beta_1,
    beta_2 = beta_2,
    sigma = sigma,
    errors = errors
  )
}


# internal ----------------------------------------------------------------

run_two_group_dense <- function(y, X1, X2, nnt) {
  # extract parameters from inputs
  n <- length(y)
  k1 <- ncol(X1)
  k2 <- ncol(X2)

  nn <- 80
  dsum <- 0.0
  dsums <- as.double(rep(-7, k1+k2+3))
  stds <- as.double(rep(-7, k1+k2+3))
  X <- cbind(X1, X2)
  fit <- .Fortran("dense_eval",as.integer(nnt), as.integer(nn), as.integer(n), as.integer(k1),
           as.integer(k2), X, y, means=dsums, dsum, stds=stds)

  out <- list(means=fit$means, stds=fit$stds)
  # TODO: get rid of data generation inside this function
}
