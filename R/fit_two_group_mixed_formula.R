#' Formula interface to `fit_two_group_mixed()`
#'
#' This function provides a formula interface to the model fit by
#' [fit_two_group_mixed()].
#'
#' @export
#' @param formula An lme4-style model formula. Currently one a single varying
#'   intercept or varying slope is supported in the "random effects" part of the
#'   model. To bypass this restriction use [fit_two_group_mixed()] directly.
#' @param data A data frame containing the variables used in the model.
#' @param ... Arguments passed to [fit_two_group_mixed()], except for `y`, `X1`,
#'   and `X2`, which are generated automatically from `data`, `fixed_formula`,
#'   and `varying_intercept`.
#' @return See [fit_two_group_mixed()].
#'
#' @examples
#' fit <- fit_two_group_mixed_formula(mpg ~ wt + as.factor(gear) + (1|cyl), mtcars)
#' fit$beta1
#' fit$beta2
#'
#' fit2 <- fit_two_group_mixed_formula_2(mtcars, mpg ~ wt + as.factor(gear), "cyl")
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

# Use lme4 to convert formula to `y`, `X1`, `X2`, trying to catch unsupported specifications.
# Currently on a single varying intercept or varying slope is supported
parse_model_formula <- function(formula, data) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Please install the lme4 package.", call. = FALSE)
  }
  lf <- lme4::lFormula(formula, data)
  bars <- lme4::findbars(lf$formula)
  if (length(bars) > 1) {
    stop("Only one varying term is currently supported.", call. = FALSE)
  }
  bars <- bars[[1]]
  pre_bar <- as.character(bars)[2] # pull out terms on lhs of bar
  if (pre_bar != "1" && !grepl("0 +", pre_bar, fixed = TRUE)) {
    # only allow terms like (1 | g) and (0 + x | g)
    stop("Currently only a varying intercept or varying slope is supported, but not both.",
         call. = FALSE)
  }
  list(
    y = stats::model.response(lf$fr),
    X1 = t(as.matrix(lf$reTrms$Zt)),
    X2 = lf$X
  )
}
