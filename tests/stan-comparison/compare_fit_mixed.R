library(cmdstanr)
library(fastNoNo)

# Generate data and parameter values --------------------------------------
set.seed(1)
n <- 20
k1 <- 10
k2 <- 10
sd_beta2 <- sample.int(0.5, k2, replace = TRUE)
sd_sigma_y <- 1.1
sd_sigma1 <- 1.2
X1 <- rnorm(n * k1, mean = 0, sd = 1)
X1 <- matrix(data = X1, nrow = n, ncol = k1)
X2 <- rnorm(n * k2, mean = 0, sd = 1)
X2 <- matrix(data = X2, nrow = n, ncol = k2)
X <- cbind(X1, X2)
beta1 <- rnorm(k1, mean = 0, sd = sd_sigma1)
beta2 <- rnorm(k2, mean = 0, sd = sd_beta2)
beta <- c(beta1, beta2)
y <- array(X %*% beta + rnorm(n, mean = 0, sd = sd_sigma_y))

# Stan posterior means and sds --------------------------------------------
data_list <- list(
  n = n, k1 = k1, k2 = k2,
  y = y, X1 = X1, X2 = X2,
  sd_beta2 = sd_beta2, sd_sigma_y = sd_sigma_y, sd_sigma1 = sd_sigma1
)
mod <- cmdstan_model("~/posterior_calib/fit_mixed.stan")
fit_stan <- mod$sample(data = data_list, seed = 1, parallel_chains = 4, iter_sampling = 1e3)
posterior_summary <- fit_stan$summary()[-1, ] # drop lp__
stan_estimates <- as.data.frame(posterior_summary[, c("mean", "sd")])
rownames(stan_estimates) <- posterior_summary$variable

# fastNoNo posterior means and sds ----------------------------------------
fit_fastnono <- fit_mixed(y, X1, X2, sd_beta2 = sd_beta2, sd_sigma_y = sd_sigma_y, sd_sigma1 = sd_sigma1, nnt = 20)
fastnono_estimates <- rbind(fit_fastnono$sigma, fit_fastnono$beta1, fit_fastnono$beta2)
rownames(fastnono_estimates) <- rownames(stan_estimates)

# differences between Stan and fastNoNo -----------------------------------
diffs <- (stan_estimates - fastnono_estimates) / fastnono_estimates
diffs
