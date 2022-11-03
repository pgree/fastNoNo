#' Formula interface to `fit_mixed()`
#'
#' This function provides a limited `lmer`-style formula interface to the model
#' fit by [fit_mixed()]. Currently only a small subset of possible formula
#' specifications are supported.
#'
#' @export
#' @param formula (formula) An `lmer`-style model formula. Currently only a
#'   single varying intercept or varying slope is supported in the "random
#'   effects" part of the model formula. To bypass this restriction use
#'   [fit_mixed()] directly.
#' @param data (data frame) The data containing the variables used in the model.
#' @param ... Currently for internal use only.
#' @inheritParams fit_mixed
#' @inherit fit_mixed return
#' @inherit fit_mixed references
#'
#' @seealso [fit_mixed()]
#' @examples
#' fit <- fit_mixed_formula(
#'   formula = mpg ~ wt + as.factor(gear) + (1|cyl),
#'   data = mtcars,
#'   sd_beta2 = 10,
#'   sd_sigma_y = 10,
#'   sd_sigma1 = 5
#' )
#' fit$beta1
#' fit$beta2
#' fit$sigma
#'
#' # check accuracy of estimates
#' fit$errors
#'
#' # refitting using more quadrature nodes improves accuracy
#' fit <- fit_mixed_formula(
#'   formula = mpg ~ wt + as.factor(gear) + (1 | cyl),
#'   data = mtcars,
#'   sd_beta2 = 10,
#'   sd_sigma_y = 10,
#'   sd_sigma1 = 5,
#'   nnt = 30
#' )
#' fit$errors
#'
fit_mixed_formula <- function(formula, data, ...,
                              sd_sigma_y, sd_sigma1, sd_beta2,
                              nnt = 10) {
  stopifnot(
    is.data.frame(data),
    !anyNA(data)
  )
  model_data <- parse_model_formula(formula, data, on_failed_check = "error")
  out <- fit_mixed(
    y = model_data$y,
    X1 = model_data$X1,
    X2 = model_data$X2,
    sd_sigma_y = sd_sigma_y,
    sd_sigma1 = sd_sigma1,
    sd_beta2 = sd_beta2,
    nnt = nnt
  )
  if (isTRUE(list(...)$debug)) { # temporary to help with debugging
    out$debug <- list(X1 = model_data$X1, X2 = model_data$X2)
  }
  out
}


# internal ----------------------------------------------------------------

#' Parse model formula and return data that can be passed to fastNoNo
#'
#' Uses lme4 to convert formula to `y`, `X1`, `X2`, trying to catch unsupported
#' specifications. Currently only a single varying intercept or varying slope is
#' supported.
#'
#' @noRd
#' @param formula,data The model formula and data frame.
#' @param on_failed_check What to do if the formula contains unsupported
#'   terms. Can be `"error"`, `"warning"`, or `"ignore"`.
#' @return A list with vector `y` and matrices `X1` and `X2`.
#'
parse_model_formula <- function(formula, data, on_failed_check = c("error", "warning", "ignore")) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Please install the lme4 package.", call. = FALSE)
  }
  lf <- lme4::lFormula(
    formula = formula,
    data = data,
    control = lme4::lmerControl(
      check.nlev.gtreq.5 = "ignore",
      check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nlev = "ignore",
      check.nobs.vs.nRE = "ignore",
      check.scaleX = "ignore"
    )
  )
  check_formula_unsupported_terms(
    formula = lf$formula,
    on_failed_check = match.arg(on_failed_check)
  )
  list(
    y = stats::model.response(lf$fr),
    X1 = make_X1_matrix(lf$reTrms),
    X2 = lf$X
  )
}

#' Make the `X1` matrix to pass to fastNoNo
#' @noRd
#' @param reTrms The `reTrms` object returned by `lme4::lFormula()`.
#' @return The `X1` matrix to pass to fastNoNo.
#'
make_X1_matrix <- function(reTrms) {
  # for now we have to convert Zt from sparse to dense in order to pass a matrix
  # (and not a dgCMatrix object) to fastNoNo, but maybe we can preserve the
  # sparse matrix in the future
  X1_t <- as.matrix(reTrms$Zt)
  X1 <- t(X1_t)
  colnames(X1) <- make_X1_colnames(reTrms)
  X1
}

#' Generate informative column names for X1 matrix
#'
#' The Z matrix in lme4 (our X1 matrix) only uses levels of the grouping
#' variable as column names, so the names will repeat if there are a
#' varying intercept and slope by the same grouping variable (e.g. `(1+x|g)`)
#' or if there are multiple varying intercepts and the grouping variables happen
#' to have overlapping level names (e.g. `(1|g) + (1|f)` and both `g` and `f`
#' share level names). In these cases `Z` will have repetitive column names,
#' which prevents using those names as parameter names when returning fastNoNo
#' estimates to the user. This function creates informative names that don't
#' repeat.
#'
#' @noRd
#'
#' @note This function isn't really needed given that we currently only allow 1
#'   varying term, but it will be useful if we eventually allow more complicated
#'   model formulas.
#'
#' @param reTrms The `reTrms` object returned by `lme4::lFormula()`.
#'
#' @return Character vector of column names to use for `X1`. Names are of the
#'   form `(Intercept)_g:level` for varying intercepts by level of grouping
#'   variable `g`, and `x_g:level` for varying slopes on a variable `x` by level
#'   of grouping variable `g`. For example, with the `mtcars` data and a model
#'   formula `~ (1|cyl) + (1+disp|gear)` we get names like `(Intercept)_cyl:4`,
#'   `(Intercept)_gear:3`, `disp_gear:3`, etc.
#'
make_X1_colnames <- function(reTrms) {
  X1_colnames <- c()
  for (i in seq_along(reTrms$Ztlist)) {
    nms <- reTrms$Ztlist[[i]]@Dimnames[[1]]
    nms <- paste0(names(reTrms$cnms)[i], ":", nms)
    nms <- paste0(rep(reTrms$cnms[[i]], times = length(unique(nms))), "_", nms)
    X1_colnames <- c(X1_colnames, nms)
  }
  X1_colnames
}

#' Check if model formula has unsupported terms
#'
#' Currently formulas are restricted to a single varying intercept or slope,
#' e.g. `(1|g)` or `(0 + x|g)`.
#'
#' @noRd
#' @param formula lmer-style model formula.
#' @param action The action to take if the check fails (error, warning, or ignore).
#' @return Returns `TRUE` invisibly unless an error is thrown.
#'
check_formula_unsupported_terms <- function(formula, on_failed_check = c("error", "warning", "ignore")) {
  action <- match.arg(on_failed_check)
  if (action == "ignore") {
    return(invisible(TRUE))
  }

  # check for multiple expressions with "|"
  bars <- lme4::findbars(formula)
  if (length(bars) > 1) {
    switch(
      action,
      warning = warn_multiple_varying_terms(),
      error = stop("Only one varying term is currently supported.", call. = FALSE)
    )
  }

  # check that all terms are of the form (1 | g) or (0 + x | g)
  # (currently bars is always length one but this may change)
  for (j in seq_along(bars)) {
    bars_j <- bars[[j]]
    pre_bar <- as.character(bars_j)[2] # pull out terms on lhs of bar
    if (pre_bar != "1" && !grepl("0 +", pre_bar, fixed = TRUE)) {
      switch(
        action,
        warning = warn_multiple_varying_terms(),
        error = stop(
          "Currently only terms (1 | g) and (0 + x | g) are supported. ",
          "(1 + x|g) is not yet implemented.",
          call. = FALSE
        )
      )
    }
  }
  invisible(TRUE)
}

warn_multiple_varying_terms <- function() {
  warning(
    "Multiple varying terms with partial pooling detected: \n",
    "Currently, we advise using this function with partial pooling for only a ",
    "single varying intercept or slope because this implementation of the fastNoNo ",
    "algorithm uses same standard deviation parameter (sigma1) for all varying terms.",
    call. = FALSE
  )
}
