
# test that formula interface returns same results as original ------------

X2 <- model.matrix(~ wt + disp + factor(gear), data = mtcars)
X1 <- model.matrix(~ 0 + as.factor(cyl), data = mtcars)
fit <- fit_mixed(mtcars$mpg, X1, X2, nnt = 100, sd_beta2 = 10)
fit_formula <- fit_mixed_formula(mpg ~ wt + disp + factor(gear) + (1|cyl), mtcars, nnt = 100, sd_beta2 = 10)

# for beta1 avoid checking rownames (they will currently differ for beta1)
expect_equal(fit_formula$beta1, fit$beta1, check.attributes = FALSE)
expect_equal(fit_formula$beta2, fit$beta2)
expect_equal(fit_formula$sigma, fit$sigma)
expect_equal(fit_formula$errors, fit$errors, check.attributes = FALSE)
expect_equal(fit_formula$cov, fit$cov, check.attributes = FALSE)


# test that errors thrown for unsupported formula terms -------------------

# valid formulas
expect_silent(
  fit_mixed_formula(mpg ~ (1 | cyl), data = mtcars)
)
expect_silent(
  fit_mixed_formula(mpg ~ (0 + wt | cyl), data = mtcars)
)

# currently invalid formulas
expect_error(
  fit_mixed_formula(mpg ~ (1 + wt | cyl), data = mtcars),
  "Currently only terms (1 | g) and (0 + x | g) are supported. (1 + x|g) is not yet implemented.",
  fixed = TRUE
)
expect_error(
  fit_mixed_formula(mpg ~ (1 | cyl/gear), data = mtcars),
  "Only one varying term is currently supported."
)
expect_error(
  fit_mixed_formula(mpg ~ (1|cyl) + (1|gear), data = mtcars),
  "Only one varying term is currently supported."
)
expect_error(
  fit_mixed_formula(mpg ~ (1 + wt|cyl) + (1|gear), data = mtcars),
  "Only one varying term is currently supported."
)

# on_failed_check = "ignore" allows formulas with unsupported terms when using parse_model_formula()
expect_silent(
  fastNoNo:::parse_model_formula(mpg ~ (1 + wt|cyl) + (1|gear), data = mtcars, on_failed_check = "ignore"),
)
expect_error(
  fastNoNo:::parse_model_formula(mpg ~ (1 + wt|cyl) + (1|gear), data = mtcars, on_failed_check = "error"),
  "Only one varying term is currently supported."
)


