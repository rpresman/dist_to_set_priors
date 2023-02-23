
# Load packages
library(Rfast)
library(rstan)
library(ggplot2)
library(forecast)
library(gridExtra)
library(bayesplot)
library(knitr)
library(kableExtra)
library(scatterplot3d)
library(lattice)
library(GGally)
library(xtable)

# Initialize parameters
F_vec <- c(1/sqrt(3), 1/sqrt(3), 1/sqrt(3))
set.seed(0)
stan_dat <- list(N = 200, F = F_vec,
                 lambda1 = 10e-5, sigma2 = 0.1, m = 5)

# Fit level set relaxation prior (rVMF)
# .stan files (here and below) adapted from:
# Leo Duan. BayesCoRe. https://github.com/leoduan/BayesCoRe. Accessed: Sept. 7, 2022.
set.seed(0)
fit_stan <- stan("relaxed_vmf_sampling.stan", data = stan_dat, chains = 2, refresh = 0)


# Fit distance-to-set relaxation prior (rVMF)
set.seed(0)
fit_stan_sqd <- stan("relaxed_sqd_vmf_sampling.stan", data = stan_dat, chains = 2, refresh = 0)

# Simulate from theoretical vMF
set.seed(0)
x_vmf <- rvmf(200, F_vec, 5)

# Extract and format output; thin by a factor of 10 (for visualization)
colnames(x_vmf) <- c("x", "y", "z")
theta_t <- as.matrix(fit_stan)[seq(10, 2000, by = 10), 1:3]
colnames(theta_t) <- c("x", "y", "z")

theta_t_sqd <- as.matrix(fit_stan_sqd)[seq(10, 2000, by = 10), 1:3]
colnames(theta_t_sqd) <- c("x", "y", "z")

x_vmf <- as.data.frame(x_vmf)
theta_t <- as.data.frame(theta_t)
theta_t_sqd <- as.data.frame(theta_t_sqd)

# Plots
p1 <- ggpairs(x_vmf, upper = "blank", diag = "blank")
p1$plots <- p1$plots[c(4:5, 7:8)]
p1$xAxisLabels <- p1$xAxisLabels[1:2]
p1$yAxisLabels <- p1$yAxisLabels[2:3]
p1

p2 <- ggpairs(theta_t, upper = "blank", diag = "blank")
p2$plots <- p2$plots[c(4:5, 7:8)]
p2$xAxisLabels <- p2$xAxisLabels[1:2]
p2$yAxisLabels <- p2$yAxisLabels[2:3]
p2

p3 <- ggpairs(theta_t_sqd, upper = "blank", diag = "blank")
p3$plots <- p3$plots[c(4:5, 7:8)]
p3$xAxisLabels <- p3$xAxisLabels[1:2]
p3$yAxisLabels <- p3$yAxisLabels[2:3]
p3

# Simulation Study
set.seed(0)
lambda_vec <- c(1e-3, 1e-4, 1e-5, 1e-6)

time_vec <- numeric(length(lambda_vec))
time_vec_sqd <- numeric(length(lambda_vec))

acc_vec <- numeric(length(lambda_vec))
acc_vec_sqd <- numeric(length(lambda_vec))

x_dat <- data.frame()

i = 1

for (lambda in lambda_vec) {
  # Set up stan parameters
  stan_dat <- list(N = 200, F = F_vec,
                   lambda1 = lambda, sigma2 = 0.1, m = 5)
  
  # Run both constraint relaxation approaches
  fit_stan <- stan("relaxed_vmf_sampling.stan",
                   data = stan_dat,chains = 2, refresh = 0)
  fit_stan_sqd <- stan("relaxed_sqd_vmf_sampling.stan",
                       data =stan_dat, chains = 2, refresh = 0)
  
  # Extract summary statistics for each approach
  fit_stan_sum <- summary(fit_stan)$summary
  fit_stan_sqd_sum <- summary(fit_stan_sqd)$summary
  
  x <- fit_stan_sum[1:3, c("mean", "2.5%", "97.5%", "n_eff")]
  x <- data.frame(c("x","y","z"), x)
  rownames(x) <- NULL
  x <- data.frame(lambda = c(1/lambda, "", ""), x)
  colnames(x) <- c("rho", "theta", "Mean",
                   "2.5%", "97.5%", "ESS")
  
  y <- fit_stan_sqd_sum[1:3, c("mean", "2.5%", "97.5%", "n_eff")]
  rownames(y) <- NULL
  colnames(y) <- c("Mean",
                   "2.5%", "97.5%", "ESS")
  
  x_dat <- rbind(x_dat, cbind(x,y))
  
  # Extract performance measures for each approach
  time_vec[i] <- mean(get_elapsed_time(fit_stan)[,2])
  time_vec_sqd[i] <- mean(get_elapsed_time(fit_stan_sqd)[,2])
  
  sampler_params <- get_sampler_params(fit_stan,
                                       inc_warmup = FALSE)
  sampler_params_sqd <- get_sampler_params(fit_stan_sqd,
                                           inc_warmup = FALSE)
  acc_vec[i] <- mean(sapply(sampler_params,
                            function(x) mean(x[, "accept_stat__"])))
  acc_vec_sqd[i] <- mean(sapply(sampler_params_sqd,
                                function(x) mean(x[, "accept_stat__"])))
  
  i = i + 1
}

# Comparison of sampling performance between level set relaxation and distance-to-set relaxation priors
x_dat %>% 
  kable(digits = 3, booktabs = TRUE, align = "c",
        col.names = c("$\\rho$", "Coordinate", "Mean", "2.5%", "97.5%", "ESS",
                      "Mean", "2.5%", "97.5%", "ESS"),
        escape = FALSE) %>% 
  kable_classic() %>% 
  add_header_above(c(" " = 2, "Level Set Relaxation" = 4, "Distance-to-Set" = 4))

# Comparison of acceptance rates between level set relaxation and distance-to-set relaxation priors
acc_dat <- data.frame(1/lambda_vec,
                      acc_vec,
                      acc_vec_sqd)
acc_dat$X1.lambda_vec <- as.character(acc_dat$X1.lambda_vec)

acc_dat %>% 
  kable(col.names = c("$\\rho$", "Level Set Relaxation", "Distance-to-Set"),
        digits = 3, booktabs = TRUE, align = "c", escape = F) %>% 
  kable_classic()