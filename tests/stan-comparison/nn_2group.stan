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
  vector[k1] alpha1;
  vector[k2] alpha2;

}
model {
  y ~ normal(X1*alpha1 + X2*alpha2, sigma_y);
  alpha1 ~ normal(0, sigma1);
  alpha2 ~ normal(0, sigma2);
  sigma_y ~ normal(0, sigs[1]);
  sigma1 ~ normal(0, sigs[2]);
  sigma2 ~ normal(0, sigs[3]);
}
