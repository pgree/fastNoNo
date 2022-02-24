data {
  int n;
  int k1;
  int k2;
  vector[n] y;
  matrix[n,k1] X1;
  matrix[n,k2] X2;
  vector[k2] ss;
  real sdy;
  real sd1;
}
parameters {
  real<lower=0> sigma_y;
  real<lower=0> sigma_beta1;
  vector[k1] beta1;
  vector[k2] beta2;
}
model {
  y ~ normal(X1*beta1 + X2*beta2, sigma_y);
  beta1 ~ normal(0, sigma_beta1);
  for (i in 1:k2) {
    beta2[i] ~ normal(0, ss[i]);
  }
  sigma_y ~ normal(0, sdy);
  sigma_beta1 ~ normal(0, sd1);
}
