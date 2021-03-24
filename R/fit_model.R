#' Fit a two-group normal-normal model
#'
#' @export
#' @param X The matrix of data of size n x (k1+k2)
#' @param k1 size of first group
#' @param k2 size of second group
#' @param y outcome
#' @return fit a list containing posterior means, standard devations,
#' and accuracy
#'
fit_model <- function(X, k1, k2, y) {
  # extract number of observations
  n <- dim(X)[1]

  # write data to temp file to be read by fortran
  tmp_dir <- tempdir()
  filename <- file.path(tmp_dir, "params.dat")
  write_data(n, k1, k2, X, y, filename)

  # run fortran
  run_fortran(tmp_dir)

  # read results
  fit <- read_output(tmp_dir)

  return(fit)
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


run_fortran <- function(dir) {
  # locations of executables
  compiled_exe <- system.file("two_group_dense_compiled", package = "fastNN")
  fortran_exe <- system.file("two_group_dense", package = "fastNN")
  f_obj_file <- system.file("two_group_dense.o", package="fastNN")
  inst_dir <- system.file("", package = "fastNN")

  # if the object file already exists, don't recompile
  if (file.exists(f_obj_file)) {
    processx::run(compiled_exe, dir, wd=inst_dir, error_on_status = FALSE)
  } else {
    processx::run(fortran_exe, dir, wd=inst_dir, error_on_status = FALSE)
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
  k <- dim(df_out)[1]-3
  df_out$cols <- c(1:k, 'sig1', 'sig2', 'sig3')

  # add columns names
  colnames(df_out)[1] <- "mean"
  colnames(df_out)[2] <- "std"

  return(c(df_out, error=derr))
}
