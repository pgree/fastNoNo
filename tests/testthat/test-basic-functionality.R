# for now just tests that it runs without erroring,
# not whether the answer is correct

test_that("fit_two_group runs", {
  n <- 10000
  k_1 <- 50
  k_2 <- 60

  sigma_y <- 1
  sigma_1 <- 0.5
  sigma_2 <- 2
  beta_1 <- rnorm(k_1, 0, sigma_1)
  beta_2 <- rnorm(k_2, 0, sigma_2)

  X_1 <- matrix(rnorm(n * k_1, 2, 3), ncol = k_1)
  X_2 <- matrix(rnorm(n * k_2, -1, 5), ncol = k_2)
  y <- rnorm(n, X_1 %*% beta_1 + X_2 %*% beta_2, sigma_y)

  # hyperpriors on scale parameters
  sd_1 <- 1
  sd_2 <- 2
  sd_y <- 2

  # Fit model
  fit <- fit_two_group(y, X_1, X_2, sd_y, sd_1, sd_2)
  expect_named(fit, c("beta_1", "beta_2", "sigma", "errorsXXXX"))
})

test_that("fit_two_group_mixed runs", {
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

  # Fit model
  fit <- fit_two_group_mixed(y, X_1, X_2, ss = rep(1, k_2), sd_y = 1, sd_1 = 1, nnt = 20)
  expect_named(fit, c("beta_1", "beta_2", "sigma", "cov", "errors"))
})
