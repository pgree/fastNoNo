#' Two-group normal-normal model fit
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
#' Currently the algorithm is implemented in Fortran and using the package
#' requires the ability to compile Fortran code. Eventually the Fortran code may
#' come pre-compiled.
#'
#' @export
#' @param X1 Data matrix corresponding to group 1.
#' @param X2 Data matrix corresponding to group 2.
#' @param y Outcome vector.
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
#' k1 <- 50
#' k2 <- 60
#' n <- 10000
#'
#' y <- rnorm(n)
#' X1 <- matrix(rnorm(n * k1, 2, 3), ncol = k1)
#' X2 <- matrix(rnorm(n * k2, -1, 5), ncol = k2)
#'
#' fit <- fit_two_group_dense(X1, X2, y)
#' str(fit)
#' }
#'
fit_two_group_dense <- function(X1, X2, y) {
  # check matrices are same length
  stopifnot(nrow(X1) == nrow(X2), length(y) == nrow(X1))

  # extract number of observations
  n <- nrow(X1)
  k1 <- ncol(X1)
  k2 <- ncol(X2)
  X <- cbind(X1, X2)

  # write data to temp file to be read by fortran
  tmp_dir <- tempdir()
  filename <- file.path(tmp_dir, "params.dat")
  write_data(n, k1, k2, X, y, filename)

  # run fortran
  run_fortran(tmp_dir)

  # read results
  read_output(tmp_dir, k1, k2)
}

write_data <- function(n, k1, k2, a, y, filename) {
  # make sure arguments are of correct type
  stopifnot(is.matrix(a))
  stopifnot(n > (k1+k2))

  # write data to file in tempdir()
  nj <- sprintf("%012d", c(n, k1, k2))
  write.table(nj, file=filename, row.names = FALSE, col.names = FALSE, sep=',', quote=FALSE)
  A2 <- sprintf("%012.6f", a)
  write.table(A2, file=filename, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
  y2 <- sprintf("%012.6f", y)
  write.table(y2, file=filename, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
}


run_fortran <- function(compile_dir) {
  # locations of executables
  create_fortran_exe <- system.file("two_group_dense", package = "fastNoNo")
  fortran_exe <- file.path(compile_dir, "int2")
  inst_dir <- system.file("", package = "fastNoNo")

  # if the executable already exists, don't recompile
  # the second argument in the command below provides the location of the
  # params.dat file that is read by the fortran
  if (!file.exists(fortran_exe)) {
    message("Compiling fortran code. This may take a moment...")
    processx::run(create_fortran_exe, compile_dir, wd=inst_dir, error_on_status = FALSE)
  }

  # run exectuable
  processx::run(fortran_exe, compile_dir, wd=compile_dir, error_on_status = FALSE)

}

read_output <- function(dir, k1, k2) {
  # read the posterior means and stds written to a text file in tmp_dir
  # by fortran
  filename <- file.path(dir, "exps.dat")

  # read results from fortran
  df_out <- read.csv(file=filename, header=FALSE)

  # get error
  err <- df_out$V1[1]

  # remove row with error
  df_out <- df_out[2:nrow(df_out), , drop=FALSE]

  # add columns names
  colnames(df_out)[1] <- "mean"
  colnames(df_out)[2] <- "sd"

  df_beta_1 <- df_out[1:k1, ]
  df_beta_2 <- df_out[(k1+1):(k1+k2), ]
  df_sigma <- df_out[(nrow(df_out)-2):nrow(df_out), ]

  # add row names
  rownames(df_beta_1) <- paste0("beta_1_", 1:k1)
  rownames(df_beta_2) <- paste0("beta_2_", 1:k2)
  rownames(df_sigma) <- c("sigma_y", "sigma_1", "sigma_2")

  list(
    beta_1 = df_beta_1,
    beta_2 = df_beta_2,
    sigma = df_sigma,
    error = err
  )
}
