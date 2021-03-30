#' Two-group normal-normal model fit
#'
#' This function fits the Bayesian two-group linear regression model
#' \deqn{y ~ normal(X_1 \beta_1 + X_2 \beta_2, \sigma_y)}
#' \deqn{\beta_1 ~ normal(0, \sigma_1)}
#' \deqn{\beta_2 ~ normal(0, \sigma_2)}
#' \deqn{\sigma_1 ~ normal+(0, 1)}
#' \deqn{\sigma_2 ~ normal+(0, 1)}
#' \deqn{\sigma_y ~ normal+(0, 1)}
#'
#' @export
#' @param X1 Data matrix corresponding to group 1.
#' @param X2 Data matrix corresponding to group 2.
#' @param y Outcome vector.
#' @return A list containing posterior means, standard deviations,
#' and the approximate accuracy of the standard deviations and means. The `mean`
#' and `std` vectors in the list are of length k+3. The first k entries
#' correspond to posterior means and standard deviations of the k parameters
#' \eqn{\beta_1, \beta_2} and the final three correspond to the scale parameters
#' \eqn{\sigma_1, \sigma_2, \sigma_y}.
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
  read_output(tmp_dir)
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


run_fortran <- function(tmp_dir) {
  # locations of executables
  fortran_exe <- system.file("two_group_dense", package = "fastNoNo")
  fortran_exe2 <- file.path(tmp_dir, "int2")
  inst_dir <- system.file("", package = "fastNoNo")

  # if the executable already exists, don't recompile
  # the second argument in the command below provides the location of the
  # params.dat file that is read by the fortran
  if (file.exists(fortran_exe2)) {
    processx::run(fortran_exe2, tmp_dir, wd=tmp_dir, error_on_status = FALSE)
  } else {
    message("Compiling fortran code. This may take a moment...")
    processx::run(fortran_exe, tmp_dir, wd=inst_dir, error_on_status = FALSE)
  }

}

read_output <- function(dir) {
  # read the posterior means and stds written to a text file in tmp_dir
  # by fortran
  filename <- file.path(dir, "exps.dat")

  # read results from fortran
  df_out <- read.csv(file=filename, header=FALSE)

  # get error
  derr <- df_out$V1[1]

  # remove row with error
  df_out <- df_out[2:nrow(df_out), , drop=FALSE]
  ###k <- dim(df_out)[1]-3
  ###df_out$cols <- c(1:k, 'sig1', 'sig2', 'sig3')

  # add columns names
  colnames(df_out)[1] <- "mean"
  colnames(df_out)[2] <- "std"

  return(c(df_out, error=derr))
}
