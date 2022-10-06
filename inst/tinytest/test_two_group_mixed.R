set.seed(123)

n <- 1000
k_1 <- 50
k_2 <- 60
sigma_y <- 1
sigma_1 <- 0.5
beta_1 <- rnorm(k_1, 0, sigma_1)
beta_2 <- rnorm(k_2, 0, 1)
X_1 <- matrix(rnorm(n * k_1, 2, 3), ncol = k_1)
X_2 <- matrix(rnorm(n * k_2, -1, 5), ncol = k_2)
y <- rnorm(n, X_1 %*% beta_1 + X_2 %*% beta_2, sigma_y)


# test input checking -----------------------------------------------------

# certain values not allowed
y_NA <- y
X_1_NA <- X_1
X_2_NA <- X_2
y_NA[1] <- NA
X_1_NA[1,1] <- NA
X_2_NA[1,1] <- NA

expect_error(
  fit_two_group_mixed(y_NA, X_1, X_2),
  "!anyNA(y) is not TRUE",
  fixed = TRUE
)
expect_error(
  fit_two_group_mixed(y, X_1_NA, X_2),
  "!anyNA(X1) is not TRUE",
  fixed = TRUE
)
expect_error(
  fit_two_group_mixed(y, X_1, X_2_NA),
  "!anyNA(X2) is not TRUE",
  fixed = TRUE
)
expect_error(
  fit_two_group_mixed(y, X_1, X_2, ss = -1),
  "all(ss > 0) is not TRUE",
  fixed = TRUE
)
expect_error(
  fit_two_group_mixed(y, X_1, X_2, sd_y = -1),
  "sd_y > 0 is not TRUE",
  fixed = TRUE
)
expect_error(
  fit_two_group_mixed(y, X_1, X_2, sd1 = -1),
  "sd1 > 0 is not TRUE",
  fixed = TRUE
)
expect_error(
  fit_two_group_mixed(y, X_1, X_2, nnt = -20),
  "nnt > 0 is not TRUE",
  fixed = TRUE
)
expect_error(
  fit_two_group_mixed(y, X_1, X_2, nnt = pi),
  "nnt == as.integer(nnt) is not TRUE",
  fixed = TRUE
)


# incorrect sizes
expect_error(
  fit_two_group_mixed(y[-1], X_1, X_2),
  "length(y) == nrow(X1) is not TRUE",
  fixed = TRUE
)
expect_error(
  fit_two_group_mixed(y, X_1[-1, ], X_2),
  "nrow(X1) == nrow(X2) is not TRUE",
  fixed = TRUE
)
expect_error(
  fit_two_group_mixed(y, X_1, X_2[-1, ]),
  "nrow(X1) == nrow(X2) is not TRUE",
  fixed = TRUE
)
expect_error(
  fit_two_group_mixed(y, X_1, X_2, ss = c(1, 2, 3)),
  "length(ss) == 1 || length(ss) == ncol(X2) is not TRUE",
  fixed = TRUE
)
expect_error(
  fit_two_group_mixed(y, X_1, X_2, sd_y = c(1, 2)),
  "length(sd_y) == 1 is not TRUE",
  fixed = TRUE
)
expect_error(
  fit_two_group_mixed(y, X_1, X_2, sd1 = c(1, 2)),
  "length(sd1) == 1 is not TRUE",
  fixed = TRUE
)
expect_error(
  fit_two_group_mixed(y, X_1, X_2, nnt = c(1, 2)),
  "length(nnt) == 1 is not TRUE",
  fixed = TRUE
)


# test contents of fit object ---------------------------------------------
fit <- fit_two_group_mixed(y, X_1, X_2)

expect_equal(names(fit), c("beta1", "beta2", "sigma", "cov", "errors", "time"))

expect_inherits(fit$beta1, "data.frame")
expect_equal(dim(fit$beta1), c(k_1, 2))
expect_equal(colnames(fit$beta1), c("mean", "sd"))
expect_equal(rownames(fit$beta1), paste0("beta1", "_", 1:k_1))

expect_inherits(fit$beta2, "data.frame")
expect_equal(dim(fit$beta2), c(k_2, 2))
expect_equal(colnames(fit$beta2), c("mean", "sd"))
expect_equal(rownames(fit$beta2), paste0("beta2", "_", 1:k_2))

expect_inherits(fit$sigma, "data.frame")
expect_equal(dim(fit$sigma), c(2, 2))
expect_equal(colnames(fit$sigma), c("mean", "sd"))
expect_equal(rownames(fit$sigma), c("sigma_y", "sigma_beta1"))

expect_inherits(fit$errors, "data.frame")
expect_equal(dim(fit$errors), c(k_1 + k_2 + 2, 2))
expect_equal(colnames(fit$errors), c("error_mean", "error_sd"))
expect_equal(rownames(fit$errors), c(rownames(fit$beta1), rownames(fit$beta2), rownames(fit$sigma)))

expect_inherits(fit$cov, "matrix")
expect_equal(dim(fit$cov), c(k_1 + k_2, k_1 + k_2))
expect_equal(rownames(fit$cov), c(rownames(fit$beta1), rownames(fit$beta2)))
expect_equal(rownames(fit$cov), colnames(fit$cov))

expect_inherits(fit$time, "numeric")
expect_equal(length(fit$time), 1)



# test algorithm is deterministic -----------------------------------------

# with the same inputs we should get the same answers if we refit
fit1 <- fit_two_group_mixed(y, X_1, X_2, ss = 1, sd_y = 1, sd1 = 1)
fit2 <- fit_two_group_mixed(y, X_1, X_2, ss = 1, sd_y = 1, sd1 = 1)

expect_equal(fit1$beta1, fit2$beta1)
expect_equal(fit1$beta2, fit2$beta2)
expect_equal(fit1$sigma, fit2$sigma)
expect_equal(fit1$errors, fit2$errors)
expect_equal(fit1$cov, fit2$cov)



# test prior arguments behave reasonably ----------------------------------

# these are just sanity checks, nothing precise

fit1 <- fit_two_group_mixed(y, X_1, X_2, sd_y = 10)
fit2 <- fit_two_group_mixed(y, X_1, X_2, sd_y = .1)
expect_true(all(fit1$sigma["sigma_y", ] > fit2$sigma["sigma_y", ]))

fit1 <- fit_two_group_mixed(y, X_1, X_2, sd1 = 10)
fit2 <- fit_two_group_mixed(y, X_1, X_2, sd1 = .1)
expect_true(all(fit1$sigma["sigma_beta1", ] > fit2$sigma["sigma_beta1", ]))



# test that results haven't changed ---------------------------------------

fit <- fit_two_group_mixed(y, X_1, X_2)
# dump("fit", file = "inst/tinytest/answers/two_group_mixed-01.R")
answer <- source("answers/two_group_mixed-01.R", local = TRUE)$value
expect_equal(fit, answer)

fit <- fit_two_group_mixed(y, X_1, X_2, ss = 2, sd_y = 2, sd1 = 2)
# dump("fit", file = "inst/tinytest/answers/two_group_mixed-02.R")
answer <- source("answers/two_group_mixed-02.R", local = TRUE)$value
expect_equal(fit, answer)

fit <- fit_two_group_mixed(y, X_1, X_2, nnt = 80)
# dump("fit", file = "inst/tinytest/answers/two_group_mixed-03.R")
answer <- source("answers/two_group_mixed-03.R", local = TRUE)$value
expect_equal(fit, answer)



# test that no NaNs in results --------------------------------------------
# https://github.com/pgree/fastNoNo/issues/13

n_runs <- 100
has_NaNs <- rep(FALSE, n_runs)
for (i in 1:n_runs) {
  fit <- fit_two_group_mixed(y, X_1, X_2)
  has_NaNs[i] <- any(sapply(fit, anyNA))
}
expect_true(sum(has_NaNs) == 0)

