data {
  int n;
  int k1;
  int k2;
  vector[n] y;
  matrix[n,k1] X1;
  matrix[n,k2] X2;
  vector[3] sigs;
}
parameters {
  real<lower=0> sigma_y;
  real<lower=0> sigma1;
  real<lower=0> sigma2;
  vector[k1] beta1;
  vector[k2] beta2;

}
model {
  y ~ normal(X1*beta1 + X2*beta2, sigma_y);
  beta1 ~ normal(0, sigma1);
  beta2 ~ normal(0, sigma2);
  sigma_y ~ normal(0, sigs[1]);
  sigma1 ~ normal(0, sigs[2]);
  sigma2 ~ normal(0, sigs[3]);
}
