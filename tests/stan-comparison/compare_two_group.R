library(cmdstanr)
library(fastNoNo)

# Generate data and parameter values --------------------------------------
n <- 200
k_1 <- 3
k_2 <- 3

sigma_y <- 1
sigma_1 <- 0.5
sigma_2 <- 2
beta_1 <- rnorm(k_1, 0, sigma_1)
beta_2 <- rnorm(k_2, 0, sigma_2)

X_1 <- matrix(rnorm(n * k_1, 2, 3), ncol = k_1)
X_2 <- matrix(rnorm(n * k_2, -1, 5), ncol = k_2)
y <- rnorm(n, X_1 %*% beta_1 + X_2 %*% beta_2, sigma_y)

# hyperpriors on scale parameters
std_y <- 2
std_1 <- 3
std_2 <- 10
sigs <- c(std_y, std_1, std_2)


# Stan posterior means and sds --------------------------------------------
mod <- cmdstan_model("tests/stan-comparison/nn_two_group.stan")
data_list <- list(n=n, k1=k_1, k2=k_2, y=y, sigs=sigs, X1=X_1, X2=X_2)
fit_stan <- mod$sample(data = data_list, seed = 1, parallel_chains = 4)
posterior_summary <- fit_stan$summary()[-1, ] # drop lp__
stan_estimates <- as.data.frame(posterior_summary[, c("mean", "sd")])
rownames(stan_estimates) <- posterior_summary$variable

# fit_two_group NOT IMPLEMENTED IN C++ YET
# # fastNoNo posterior means and sds ----------------------------------------
# fit_fastnono <- fit_two_group(y, X_1, X_2, sigs[1], sigs[2], sigs[3])
# fastnono_estimates <- rbind(fit_fastnono$sigma, fit_fastnono$beta1, fit_fastnono$beta2)
# rownames(fastnono_estimates) <- rownames(stan_estimates)
#
# # differences between Stan and fastNoNo -----------------------------------
# (diffs <- stan_estimates - fastnono_estimates)

