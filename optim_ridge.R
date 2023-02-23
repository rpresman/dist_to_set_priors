
# Load packages
library(CVXR)
library(ggplot2)
library(MASS)
library(rstan)
library(forecast)
library(gridExtra)
library(bayesplot)
library(plotrix)
library(ggforce)

# Generate a true beta
set.seed(123)

n <- 100
p <- 2
beta <- rnorm(p, 0, 100)
beta <- 1.4 * beta / norm(beta, "2")

# Generate design matrix X
X <- matrix(rnorm(n*p), nrow=n)

# Generate y values
Y <- X %*% beta + rnorm(n, sd = 3)

# Report 2-norm of true beta
beta
norm(beta, "2")

# Set up stan data
stan_dat <- list(N = n, p = p,
                 X = X, y = Y, rho = 1e3)

# Fit unsquared model and extract trace plots/ACF plots
set.seed(123)
fit_stan <- stan("hmc2_linear.stan", data = stan_dat, chains = 2, refresh = 0)
beta_t <- as.matrix(fit_stan)[, 1:2]

p1 <- mcmc_trace(beta_t, pars = "beta[1]") +
  labs(y = expression(beta[1])) +
  ylim(-1, -0.2)
q1 <- ggAcf(beta_t[, 1]) +
  labs(title = expression("ACF of "*beta[1])) +
  ylim(0, 1)

# Fit squared model and extract trace plots/ACF plots
set.seed(123)
fit_stan_sqd <- stan("hmc1_squared.stan", data = stan_dat, chains = 2, refresh = 0)
beta_t_sqd <- as.matrix(fit_stan_sqd)[, 1:2]

p3 <- mcmc_trace(beta_t_sqd, pars = "beta[1]") +
  labs(y = expression(beta[1])) +
  ylim(-1, -0.2)

q3 <- ggAcf(beta_t_sqd[, 1]) +
  labs(title = expression("ACF of "*beta[1])) + 
  ylim(0, 1)

# Make summary plots
grid.arrange(p1, p3, q1, q3,
             nrow = 2, ncol = 2,
             top = "Sampling Performance for Distance (Left) and Squared Distance (Right) Priors")
r1 <- arrangeGrob(p1, p3, q1, q3,
                  nrow = 2, ncol = 2,
                  top = "Sampling Performance for Distance (Left) and Squared Distance (Right) Priors")

beta_t <- data.frame(beta_t)
ggplot() +
  geom_point(data = beta_t, aes(x = beta.1., y = beta.2.),
             color = "black", alpha = 0.1) +
  geom_point(data = NULL, aes(x = beta[1], y = beta[2]),
             color = "red") +
  geom_circle( aes(x0 = 0, y0 = 0, r = 1),
               color = "blue") +
  coord_fixed() +
  xlim(-2, 2) +
  ylim(-1.5, 1.5) +
  labs(x = expression(beta[1]),
       y = expression(beta[2]),
       title = "Sample Draws from Posterior") +
  theme(plot.title = element_text(size=22))