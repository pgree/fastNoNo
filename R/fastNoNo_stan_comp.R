library(tidyverse)
library(cmdstanr)
library(rstanarm)
library(rstan)
library(fastNoNo)

# generate random data
n <- 10
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
std_1 <- 1
std_2 <- 10
sigs <- c(std_y, std_1, std_2)

# run the stan model
model1 <- cmdstan_model("~/fastNoNo/R/nn_2group.stan")
niters <- 1e4

# parameters
set.seed(1)
data_list <- list(n=n, k1=k1, k2=k2, y=y, sigs=sigs, X1=X_1, X2=X_2)

# sample
fit_stan <- model1$sample(
  data = data_list,
  seed = 1,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = niters
)
means <- fit_stan$summary() %>% pull(mean)
means <- c(means[5:length(means)], means[2:4])
sds <- fit_stan$summary() %>% pull(sd)
sds <- c(sds[5:length(sds)], sds[2:4])
fit_stan$summary()

# fastNoNo
fit <- fit_two_group(y, X_1, X_2, sigs[1], sigs[2], sigs[3])
str(fit)

# differences between Stan and fastNoNo
mean_diffs <- means - c(fit$beta_1$mean, fit$beta_2$mean, fit$sigma$mean)
sd_diffs <- sds - c(fit$beta_1$sd, fit$beta_2$sd, fit$sigma$sd)
mean_diffs
sd_diffs

