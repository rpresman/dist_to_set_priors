
# Load packages
library(DirichletReg)
library(quadprog)
library(ggplot2)
library(dplyr)
library(coda)
library(kableExtra)
library(xtable)

# Set up observed contingency table data 
# Source (Lines 14-36):
# Deborshee Sen. Posterior_Projections_DS. https://github.com/deborsheesen/posterior_projections_DS. Accessed: Sept. 7, 2022.

k = c(59,48,44,43,25, 21,14,4,46,44, 54,49,48,47,64, 58,32,30,31,41)
X = matrix(k, 4,5)  #note: **NOT** byrow=T

E = matrix(0,4,20)   # Equality restriction matrix E
for (i in 1:4) {E[i,seq(i,20,4)] = 1}

IE = matrix(0,12,20) # Inequality restriction matrix E
for (i in 1:3){
  for (j in 1:4){
    m = matrix(0,4,5)
    m[i,1:j] = 1
    m[(i+1),1:j] = -1
    ind = 4*(i-1)+j
    IE[ind,] = c(m)
  }
}
A = rbind(E, -IE) # Total matrix
bvec = c(rep(1,4),rep(0,12))


Dmat = diag(length(k))
Amat = t(A)
a = 1

# Create functions for HMC implementation
# Code adapted from: Neal M. Radford, Simple Implementation of Hamiltonian Monte Carlo (2010).
# https://www.cs.toronto.edu/~radford/ham-mcmc-simple. Accessed: Sept. 7, 2022

# Function to evaluate the negative potential (=log of posterior)
U <- function(theta, y, gamma, rho){
  I <- dim(y)[1]; J <- dim(y)[2]
  U_theta <- sum(ddirichlet(theta, matrix(gamma, I, J), log = TRUE))
  for (i in 1:I) {
    U_theta <- U_theta +
      dmultinom(y[i,], prob = theta[i,], log = TRUE)
  }
  sol <- solve.QP(Dmat, c(theta), Amat, bvec, meq = 4)$solution
  sol <- matrix(sol, 4, 5)
  dist_sq <- norm(theta - sol, type = "F")^2
  U_theta <- U_theta - (rho/2) * dist_sq
  return(U_theta)
}

# Function to evaluate the gradient of the negative potential
gradient_posterior <- function(theta, y, gamma, rho){
  I <- dim(y)[1]; J <- dim(y)[2]
  gradient = matrix(NA, I, J)
  for (i in 1:I) {
    x <- as.vector(y[i,])
    pr <- as.vector(theta[i,])
    gradient[i,] <- (x + gamma - 1)/pr - (tail(x, 1) + gamma - 1)/(1-sum(head(pr, -1)))
  }
  sol <- solve.QP(Dmat, c(theta), Amat, bvec, meq = 4)$solution
  sol <- matrix(sol, 4, 5)
  dist <- theta - sol
  gradient <- gradient - rho * dist
  return(gradient)
}

# Function to perform leapfrog integrator to solve PDEs
leapfrog <- function(p, theta, epsilon, L, y, gamma, rho){
  I <- dim(y)[1]; J <- dim(y)[2]

  p <- p + epsilon* gradient_posterior(theta, y, gamma, rho)/2
  
  # Update theta and ensure it is in the simplex
  for (l in 1:L){
    theta_ = theta + epsilon * p
    if(sum(theta_[,-J] <= 0) > 0 | sum(rowSums(theta_[,-J]) > 1) > 0) {
      p = -p
      return(list("theta" = theta, "p" = p))
    }
    theta = theta_
    
    if(l != L){
      p = p + epsilon * gradient_posterior(theta, y, gamma, rho)
    }
  }
  
  p = p + epsilon* gradient_posterior(theta, y, gamma, rho)/2
  p = -p
  
  return(list("theta" = theta, "p" = p))
}

# Function to run HMC
HMC = function(y, gamma, R = 10, burn_in = 0, epsilon, L, verbose = TRUE, rho){
  I <- dim(y)[1]; J <- dim(y)[2]
  THETA = list()
  ACC = rep(0, R)
  
  # Initialize theta
  theta <- matrix(rep(rep(1,I)/J,J), I, J)
  
  # Run the leapfrog integrator
  for(r in 1:(burn_in + R)){
    if(r%%1000 == 0 & r<=burn_in & verbose == TRUE){
      print(paste("Warm-up: iteration:", r))
    } else if (r%%1000 == 0 & r>burn_in & verbose == TRUE){
      print(paste("Sampling: iteration:", r-burn_in))
    }
    
    # Sample momentum with small variance
    p = matrix(rnorm(I*J, 0, 0.2), I, J)
    
    # Run the leapfrog
    proposals =  leapfrog(p, theta, epsilon, L, y, gamma, rho)
    theta_new = proposals$theta
    
    # Ensure all distributions theta sum up to 1
    for (i in 1:I) {
      theta_new[i, J] <- 1 - sum(head(theta_new[i,], -1))
    }
    p_new = proposals$p
    
    # Reject automatically if any theta is outside the simplex
    if((sum(theta_new <= 0) > 0) | (sum(theta_new > 1) > 0)) {
      acc = 0
    }
    else {
      # Compute acceptance probability
      log_rho = U(theta, y, gamma, rho) -
        U(theta_new, y, gamma, rho) +
        norm(p, type = "F")^2/2 -
        norm(p_new, type = "F")^2/2
      r_unif = log(runif(1))
      
      # Accept-reject step
      if(r_unif < min(0, log_rho)){
        theta = theta_new
        acc = 1
      } else{
        acc = 0
      }
    }
    
    if(r > burn_in){
      THETA[[r - burn_in]] = theta
      ACC[r - burn_in] = acc
    }
  }
  return(list("posterior_samples" = THETA, "acceptance" = ACC))
}

# Run HMC Code
set.seed(0)
gamma = 0.1
epsilon = 0.001
rho = 750000
L = 5
R = 100000
burn_in = 1000
out = HMC(X, gamma, R, burn_in, epsilon, L, rho = rho)

# Thin out by a factor of 10
thin_factor <- 10
theta_vec = numeric(R/thin_factor)
t_ <- seq(thin_factor,R,by=thin_factor)
for (i in 1:(R/thin_factor)) {
  theta_vec[i] = out$posterior_samples[[t_[i]]][1,1]
}

# Plot of theta_{1,1} performance
theta_vec <- data.frame(iteration = 1:length(theta_vec), theta = theta_vec)
ggplot(theta_vec, aes(x = iteration, y = theta)) + 
    geom_line() + 
    labs(x = "Iteration", y = expression(theta[11]),
         title = expression("Traceplot of "*theta[11]))

# Compute credible intervals
df1 <- data.frame()
df2 <- data.frame()
df3 <- data.frame()
df4 <- data.frame()

out2 <- list()

thin_factor <- 10
t_ <- seq(thin_factor,R,by=thin_factor)
for (r in 1:(R/thin_factor)) {
  z <- t(apply(out$posterior_samples[[t_[r]]], MARGIN = 1, FUN = cumsum))
  out2[[r]] <- z
  df1 <- rbind(df1,
               cbind(1:5, z[1,], r))
  df2 <- rbind(df2,
               cbind(1:5, z[2,], r))
  df3 <- rbind(df3,
               cbind(1:5, z[3,], r))
  df4 <- rbind(df4,
               cbind(1:5, z[4,], r))
}

sum_stats <- apply(simplify2array(out2), 1:2, quantile, prob = c(0.025, 0.5, 0.975))
(ss1 <- sum_stats[1,,])
(ss2 <- sum_stats[2,,])
(ss3 <- sum_stats[3,,])

# Credible interval plot
df <- data.frame()
dosage <- c("Placebo", "Low Dose", "Medium Dose", "High Dose")
outcome <- c("Death", "Vegetative State", "Major Disability",
             "Minor Disability", "Good Recovery")

for (i in 1:4) {
  df <- rbind(df,
              cbind(1:5, ss1[i,], ss2[i,], ss3[i,], dosage[i]))
}
colnames(df) <- c("Outcome", "LB", "Median", "UB", "Dosage")
df$Outcome <- as.numeric(df$Outcome)
df$Dosage <- factor(df$Dosage, levels = c("Placebo", "Low Dose", "Medium Dose", "High Dose"))
df$LB <- as.numeric(df$LB)
df$Median <- as.numeric(df$Median)
df$UB <- as.numeric(df$UB)

ggplot(df, aes(color = Dosage)) +
    geom_step(aes(x = Outcome, y = LB), linetype = "dashed") +
    geom_step(aes(x = Outcome, y = UB), linetype = "dashed") +
    labs(y = "Cumulative Probabilities", title = "Credible Intervals for Contingency Table Probabilities")








