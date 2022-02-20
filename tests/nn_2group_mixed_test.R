library(tidyverse)
library(cmdstanr)
library(rstanarm)
library(rstan)
library(fastNoNo)

# cd to script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# generate data
# parameters
set.seed(1)
n <- 100
k1 <- 10
k2 <- 10
X1 <- rnorm(n*k1, mean=0, sd=1)
X1 <- matrix(data=X1, nrow=n, ncol=k1)
X2 <- rnorm(n*k2, mean=0, sd=1)
X2 <- matrix(data=X2, nrow=n, ncol=k2)
X <- cbind(X1, X2)
y <- rnorm(n, mean=0, sd=1)
ss <- sample.int(3, k2, replace=TRUE)
sdy <- 1.0
sd1 <- 2.0

file1 <- "nn_2group_mixed.stan"
mod1 <- cmdstan_model(file1)
niters <- 1e3

# parameters
set.seed(1)
data_list <- list(n=n, k1=k1, k2=k2, y=y, ss=ss, sdy=sdy, sd1=sd1, X1=X1, X2=X2)

# sample
fit1 <- mod1$sample(
  data = data_list,
  seed = 1,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = niters
)
means <- fit1$summary() %>% pull(mean)
means <- c(means[4:length(means)], means[2:3])
sds <- fit1$summary() %>% pull(sd)
sds <- c(sds[4:length(sds)], sds[2:3])
fit1$summary()

# fastNoNo
fit <- fit_two_group_mixed(y, X1, X2, ss = rep(1, k2), sd_y = sdy, sd_1 = sd1, nnt = 20)
str(fit)

# differences between Stan and fastNoNo
mean_diffs <- means - c(fit$beta_1$mean, fit$beta_2$mean, fit$sigma$mean)
sd_diffs <- sds - c(fit$beta_1$sd, fit$beta_2$sd, fit$sigma$sd)
mean_diffs
sd_diffs

