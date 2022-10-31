#' Formula interface to `fit_two_group_mixed()`
#'
#' This function provides a formula interface to the model fit by
#' [fit_two_group_mixed()].
#'
#' @export
#' @param formula An `lmer`-style model formula. Currently only a single varying
#'   intercept or varying slope is supported in the "random effects" part of the
#'   model. To bypass this restriction use [fit_two_group_mixed()] directly.
#' @param data A data frame containing the variables used in the model.
#' @param ... Arguments passed to [fit_two_group_mixed()], except for `y`, `X1`,
#'   and `X2`, which are generated automatically.
#' @return See [fit_two_group_mixed()].
#'
#' @examples
#' fit <- fit_two_group_mixed_formula(
#'   formula = mpg ~ wt + as.factor(gear) + (1|cyl),
#'   data = mtcars,
#'   ss = 10,
#'   sd_y = 10,
#'   sd1 = 5
#' )
#' fit$beta1
#' fit$beta2
#'
#' fit2 <- fit_two_group_mixed_formula_2(
#'   data = mtcars,
#'   fixed_formula = mpg ~ wt + as.factor(gear),
#'   varying_intercept = "cyl",
#'   ss = 10,
#'   sd_y = 10,
#'   sd1 = 5
#' )
#' fit2$beta1
#' fit2$beta2
#'
fit_two_group_mixed_formula <- function(formula, data, ...) {
  stopifnot(
    is.data.frame(data),
    !anyNA(data)
  )
  dots <- list(...)
  if (!is.null(dots$y) || !is.null(dots$X1) || !is.null(dots$X2)) {
    stop("'y', 'X1', and 'X2' should not be specified.", call. = FALSE)
  }
  model_data <- parse_model_formula(formula, data)
  out <- fit_two_group_mixed(model_data$y, model_data$X1, model_data$X2, ...)
  out$debug <- list(X1 = model_data$X1, X2 = model_data$X2) # temporary to help with debugging
  out
}

#' @rdname fit_two_group_mixed_formula
#' @export
#' @param fixed_formula (formula) A formula in the style of [stats::lm()] that
#'   specifies the outcome variable on the left side and the "fixed effects"
#'   terms on the right side. This formula is used to construct the `X2` matrix
#'   for [fit_two_group_mixed()].
#' @param varying_intercept (string) The name of the grouping variable by which
#'   the intercept should vary (with partial pooling). Currently only a single
#'   varying intercept is supported via this argument. This variable is used to
#'   construct the `X1` matrix for [fit_two_group_mixed()].
#'
fit_two_group_mixed_formula_2 <- function(data, fixed_formula, varying_intercept, ...) {
  stopifnot(
    is.data.frame(data),
    !anyNA(data),
    inherits(fixed_formula, "formula"),
    is.character(varying_intercept),
    length(varying_intercept) == 1
  )
  dots <- list(...)
  if (!is.null(dots$y) || !is.null(dots$X1) || !is.null(dots$X2)) {
    stop("'y', 'X1', and 'X2' should not be specified.", call. = FALSE)
  }
  if (!is.factor(data[[varying_intercept]])) {
    data[[varying_intercept]] <- as.factor(data[[varying_intercept]])
  }
  X1 <- stats::model.matrix(as.formula(paste("~ 0 + ", varying_intercept)), data = data)
  X2 <- stats::model.matrix(fixed_formula, data = data)
  y <- data[[as.character(fixed_formula[[2]])]]
  out <- fit_two_group_mixed(y, X1, X2, ...)
  out$debug <- list(X1 = X1, X2 = X2) # temporary to help with debugging
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
#' @param ... Arguments passed to `lme4::lFormula()`.
#' @param strict If `TRUE` (the default) then an error is thrown if the model
#'   formula includes unsupported terms.
#' @return A list with vector `y` and matrices `X1` and `X2`.
#'
parse_model_formula <- function(formula, data, ..., strict = TRUE) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Please install the lme4 package.", call. = FALSE)
  }
  lf <- lme4::lFormula(
    formula = formula,
    data = data,
    control = make_lmer_control(),
    ...
  )
  if (strict) {
    check_formula_unsupported_terms(lf$formula)
  }
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


#' Create model and data checking specifications for lme4 formula parser
#'
#' @noRd
#' @return A list that can be passed to the `control` argument of `lme4::lFormula()`.
#'
make_lmer_control <- function() {
  lme4::lmerControl(
    check.nlev.gtreq.5 = "ignore",
    check.nobs.vs.rankZ = "ignore",
    check.nobs.vs.nlev = "ignore",
    check.nobs.vs.nRE = "ignore",
    check.scaleX = "ignore"
  )
}


#' Error if model formula has unsupported terms
#'
#' Currently formulas are restricted to a single varying intercept or slope,
#' e.g. `(1|g)` or `(0 + x|g)`.
#'
#' @noRd
#' @return Errors if formula is not ok, otherwise returns `TRUE` invisibly.
#'
check_formula_unsupported_terms <- function(formula) {
  bars <- lme4::findbars(formula)
  if (length(bars) > 1) {
    stop("Only one varying term is currently supported.", call. = FALSE)
  }
  bars <- bars[[1]]
  pre_bar <- as.character(bars)[2] # pull out terms on lhs of bar
  if (pre_bar != "1" && !grepl("0 +", pre_bar, fixed = TRUE)) {
    # only allow terms like (1 | g) and (0 + x | g)
    stop("Currently only a single varying intercept or varying slope is supported, but not both.",
         call. = FALSE)
  }
  invisible(TRUE)
}
