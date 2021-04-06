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
#' * `error`: A scalar. The approximate accuracy of the posterior mean and
#' standard deviation estimates.
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
fit_two_group_dense <- function(y, X1, X2) {
  stopifnot(nrow(X1) == nrow(X2), length(y) == nrow(X1))

  # hard code these for now but maybe we want to expose them
  compile_dir <- tempdir()  # where to compile
  output_dir <- tempdir()   # where to write data and results

  write_data(
    n = nrow(X1),
    k1 = ncol(X1),
    k2 = ncol(X2),
    X = cbind(X1, X2),
    y = y,
    file = file.path(output_dir, "params.dat")
  )
  run_fortran(compile_dir, output_dir)
  read_output(output_dir, k1 = ncol(X1), k2 = ncol(X2))
}


# internal ----------------------------------------------------------------

write_data <- function(n, k1, k2, X, y, file) {
  stopifnot(is.matrix(X), n > (k1+k2))
  write.table(
    sprintf("%012d", c(n, k1, k2)),
    file = file,
    row.names = FALSE,
    col.names = FALSE,
    sep = ',',
    quote = FALSE
  )
  write.table(
    sprintf("%012.6f", X),
    file = file,
    append = TRUE,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  write.table(
    sprintf("%012.6f", y),
    file = file,
    append = TRUE,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}

compile_fortran <- function(compile_dir, quiet = FALSE) {
  if (!quiet) {
    message(
      "Compiling fortran code before fitting the model. " ,
      "This is only necessary the first time the function ",
      "is run in an R session."
    )
  }
  processx::run(
    command = system.file("two_group_dense", package = "fastNoNo"),
    args = compile_dir,
    wd = system.file("", package = "fastNoNo"),
    error_on_status = FALSE
  )
}

run_fortran <- function(compile_dir, output_dir) {
  fortran_exe <- file.path(compile_dir, "int2")
  if (!file.exists(fortran_exe)) {
    compile_fortran(compile_dir)
  }
  processx::run(
    command = fortran_exe,
    args = output_dir,
    wd = compile_dir,
    error_on_status = FALSE
  )

}

read_output <- function(output_dir, k1, k2) {
  df_out <- read.csv(file = file.path(output_dir, "exps.dat"), header = FALSE)

  error <- df_out$V1[1]

  estimates <- df_out[2:nrow(df_out), , drop=FALSE]
  colnames(estimates) <- c("mean", "sd")

  beta_1 <- estimates[1:k1, ]
  rownames(beta_1) <- paste0("beta_1_", 1:k1)

  beta_2 <- estimates[(k1+1):(k1+k2), ]
  rownames(beta_2) <- paste0("beta_2_", 1:k2)

  sigma <- estimates[(nrow(estimates)-2):nrow(estimates), ]
  rownames(sigma) <- c("sigma_y", "sigma_1", "sigma_2")

  list(
    beta_1 = beta_1,
    beta_2 = beta_2,
    sigma = sigma,
    error = error
  )
}
