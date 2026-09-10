## ---------------------------------------------------------------------------
## Tables 1 and 2
## Model 2: Poisson-hurdle mixture model with a bivariate normal random effect
##   Y_t = Z_t * (1 + N_t),
##     Z_t | Theta ~ Ber(logistic(Theta1)),  N_t | Theta ~ Pois(exp(Theta2)),
##     (Theta1, Theta2) ~ MVN((mu1, mu2), Sigma).
## ---------------------------------------------------------------------------
 
library(MASS)

# Sigmoid function
sigma_fun <- function(x) 1 / (1 + exp(-x))

# Conditional means and Monte Carlo standard errors
conditional_summary <- function(y1, y2) {
  y2_given_y1_eq_0 <- y2[y1 == 0]
  y2_given_y1_eq_1 <- y2[y1 == 1]
  
  if(length(y2_given_y1_eq_0) < 2 || length(y2_given_y1_eq_1) < 2) {
    stop("At least one conditioning group contains fewer than two observations.")
  }
  
  mean0 <- mean(y2_given_y1_eq_0)
  mean1 <- mean(y2_given_y1_eq_1)
  mcse0 <- sd(y2_given_y1_eq_0) / sqrt(length(y2_given_y1_eq_0))
  mcse1 <- sd(y2_given_y1_eq_1) / sqrt(length(y2_given_y1_eq_1))
  
  return(c(mean0 = mean0, mcse0 = mcse0,
           mean1 = mean1, mcse1 = mcse1,
           n0 = length(y2_given_y1_eq_0),
           n1 = length(y2_given_y1_eq_1)))
}

# Model 2 simulation
simulate_once <- function(mu1 = 0, mu2 = 0,
                          sigma1_sq = 1, sigma2_sq = 1,
                          rho = 0.5, N = 2000000) {
  
  sigma1 <- sqrt(sigma1_sq)
  sigma2 <- sqrt(sigma2_sq)
  Sigma  <- matrix(c(sigma1^2, rho*sigma1*sigma2,
                    rho*sigma1*sigma2, sigma2^2),
                  nrow = 2, byrow = TRUE)
  
  # Theta = (Theta1, Theta2) ~ MVN((mu1,mu2), Sigma)
  theta_mat <- mvrnorm(n = N, mu = c(mu1,mu2), Sigma = Sigma)
  theta1 <- theta_mat[,1]
  theta2 <- theta_mat[,2]
  
  hurdle_probability <- sigma_fun(theta1)
  poisson_mean <- exp(theta2)
  
  # Y1 = Z1 * (1 + N1)
  z1 <- rbinom(N, size = 1, prob = hurdle_probability)
  n1 <- rpois(N, lambda = poisson_mean)
  y1 <- z1 * (1 + n1)
  
  # Generate Y2 independently of Y1 conditional on the random effects
  z2 <- rbinom(N, size = 1, prob = hurdle_probability)
  n2 <- rpois(N, lambda = poisson_mean)
  y2 <- z2 * (1 + n2)
  
  return(conditional_summary(y1, y2))
}

## ---------------------------------------------------------------------------
## Table 1: sigma1^2 decreases
## ---------------------------------------------------------------------------

sigma1_sq_values <- c(5.0, 2.0, 1.0, 0.1)
results_sigma1 <- data.frame(
  sigma1_sq = sigma1_sq_values,
  mean_y2_given_y1_eq_0 = NA,
  mcse_y2_given_y1_eq_0 = NA,
  mean_y2_given_y1_eq_1 = NA,
  mcse_y2_given_y1_eq_1 = NA,
  number_y1_eq_0 = NA,
  number_y1_eq_1 = NA
)

set.seed(123) 
for(i in seq_along(sigma1_sq_values)) {
  out <- simulate_once(
    mu1 = 0,
    mu2 = 0,
    sigma1_sq = sigma1_sq_values[i],
    sigma2_sq = 1,
    rho = 0.5
  )
  results_sigma1$mean_y2_given_y1_eq_0[i] <- out["mean0"]
  results_sigma1$mcse_y2_given_y1_eq_0[i] <- out["mcse0"]
  results_sigma1$mean_y2_given_y1_eq_1[i] <- out["mean1"]
  results_sigma1$mcse_y2_given_y1_eq_1[i] <- out["mcse1"]
  results_sigma1$number_y1_eq_0[i] <- out["n0"]
  results_sigma1$number_y1_eq_1[i] <- out["n1"]
}

print(round(results_sigma1, 4))



## ---------------------------------------------------------------------------
## Table 2: sigma2^2 increases
## ---------------------------------------------------------------------------

sigma2_sq_values <- c(0.01, 0.10, 1.00, 2.00)
results_sigma2 <- data.frame(
  sigma2_sq = sigma2_sq_values,
  mean_y2_given_y1_eq_0 = NA,
  mcse_y2_given_y1_eq_0 = NA,
  mean_y2_given_y1_eq_1 = NA,
  mcse_y2_given_y1_eq_1 = NA,
  number_y1_eq_0 = NA,
  number_y1_eq_1 = NA
)

set.seed(456) 
for(i in seq_along(sigma2_sq_values)) {
  out <- simulate_once(
    mu1 = 0,
    mu2 = 0,
    sigma1_sq = 1,
    sigma2_sq = sigma2_sq_values[i],
    rho = 0.5
  )
  results_sigma2$mean_y2_given_y1_eq_0[i] <- out["mean0"]
  results_sigma2$mcse_y2_given_y1_eq_0[i] <- out["mcse0"]
  results_sigma2$mean_y2_given_y1_eq_1[i] <- out["mean1"]
  results_sigma2$mcse_y2_given_y1_eq_1[i] <- out["mcse1"]
  results_sigma2$number_y1_eq_0[i] <- out["n0"]
  results_sigma2$number_y1_eq_1[i] <- out["n1"]
}

print(round(results_sigma2, 4))



