// Adapted from: Leo Duan's GitHub BayesCoRe repo
// URL: https://github.com/leoduan/BayesCoRe/blob/master/sphere.stan

data {
  int<lower=1> N;
  vector[3] F;
  real lambda1;
  real sigma2;
  real m; //degrees of freedom in t
}

parameters {
  //vector[3] theta_vmf;
  vector[3] theta_t;
}

transformed parameters {
  vector[3] theta_t_proj;
  real<lower=0> theta_norm = 0;
  for (i in 1:3) {
    theta_norm += theta_t[i]^2;
  }
  if (theta_norm != 0)
    theta_t_proj = theta_t / sqrt(theta_norm);
}

model {
    //target +=  theta_vmf' * F /sigma2;
    
    target +=  -(m+3)/2* log(
      1+ (1+F'*F)/m/sigma2
      - 2* theta_t' * F/m/sigma2);

    //CORE: \theta  near sphere
    //target += - fabs(theta_vmf' * theta_vmf -1) / lambda1;
    target += -dot_self(theta_t - theta_t_proj) / lambda1;
}
