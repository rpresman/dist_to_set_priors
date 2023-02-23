
data {
  int<lower=0> N;
  int<lower=0> p;
  matrix[N, p] X;
  matrix[N, 1] y;
  real<lower=0> rho;
}

parameters {
  vector[p] beta;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] mu = X * beta;
  
  vector[p] beta_proj;
  real<lower=0> beta_norm;
  beta_norm = 0;
  for (i in 1:p) {
    beta_norm += beta[i]^2;
  }
  if (beta_norm >= 1)
    beta_proj = beta / sqrt(beta_norm);
  else
    beta_proj = beta;
}

model {
  sigma ~ inv_gamma(1,1);
  for (i in 1:N) {
    //target += normal_lpdf(y[i, ] | mu[i], sigma) +
    //normal_lpdf(beta | beta_proj, sigma/rho);
    target += normal_lpdf(y[i, ] | mu[i], sigma) +
    (-rho * dot_self(beta - beta_proj)) / (2 * sigma^2);
  }
}

