// Source: Leo Duan's GitHub BayesCoRe repo
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
model {
    //target +=  theta_vmf' * F /sigma2;
    
    target +=  -(m+3)/2* log(
      1+ (1+F'*F)/m/sigma2
      - 2* theta_t' * F/m/sigma2);

    //CORE: \theta  near sphere
    //target += - fabs(theta_vmf' * theta_vmf -1) / lambda1;
    target += - fabs(theta_t' * theta_t -1) / lambda1;
}
