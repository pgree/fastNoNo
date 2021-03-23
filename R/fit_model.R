#' Fit a model
#'
#' @export
#' @param data The data to use as a data frame.
#' @return A useful R object.
#'
fit_model <- function(data, ...) {
  # write data to temp file
  data_path <- write_data(data)

  # run fortran
  output_path <- run_fortran(data_path)

  # read results
  read_output(output_path)
}

write_data <- function(data) {
  # write data to file in tempdir()
  return(data_path)
}

run_fortran <- function(data_path) {
  # copy fortran file to temp dir
  dir <- tempdir()
  temp_fortran <- file.path(dir, "temp_hello.txt")
  fortran <- system.file("hello.txt", package = "fastNN")
  file.copy(fortran, temp_fortran)

  # compile in temp directory if necessary

  # execute fortran
  return(output_path)
}

read_output <- function(output_path) {
  # read in output files
  # make useful R object
  return(useful_R_object)
}
