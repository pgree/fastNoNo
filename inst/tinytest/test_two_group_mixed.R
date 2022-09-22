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
fit <- fit_two_group_mixed(y, X_1, X_2, ss = rep(1, k_2), sd_y = 1, sd1 = 1, nnt = 20)

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
expect_equal(colnames(fit$errors), c("error_means", "error_sds"))
expect_equal(rownames(fit$errors), c(rownames(fit$beta1), rownames(fit$beta2), rownames(fit$sigma)))

expect_inherits(fit$cov, "matrix")
expect_equal(dim(fit$cov), c(k_1 + k_2, k_1 + k_2))

expect_inherits(fit$time, "numeric")
expect_equal(length(fit$time), 1)

