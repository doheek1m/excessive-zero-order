## ---------------------------------------------------------------------------
## Supplementary Table S.1
## Model 2 (Poisson-hurdle mixture, bivariate normal random effect):
##   comparison of E[Y2 | Y1 = 1] and E[Y2 | Y1 = 2] as rho decreases,
##   with sigma1^2 = sigma2^2 = 0.5, mu1 = 0, mu2 = -2.
## ---------------------------------------------------------------------------
 
library(MASS)

# Sigmoid function
sigma_fun <- function(x) 1 / (1 + exp(-x))

# Conditional means and Monte Carlo standard errors
conditional_summary <- function(y1, y2) {
  y2_given_y1_eq_1 <- y2[y1 == 1]
  y2_given_y1_eq_2 <- y2[y1 == 2]
  
  mean1 <- mean(y2_given_y1_eq_1)
  mean2 <- mean(y2_given_y1_eq_2)
  mcse1 <- sd(y2_given_y1_eq_1) / sqrt(length(y2_given_y1_eq_1))
  mcse2 <- sd(y2_given_y1_eq_2) / sqrt(length(y2_given_y1_eq_2))
  
  return(c(mean1 = mean1, mcse1 = mcse1,
           mean2 = mean2, mcse2 = mcse2))
}

# Model 2 simulation
simulate_once <- function(rho, sigma1_sq = 0.5, sigma2_sq = 0.5,
                          mu1 = 0, mu2 = -2, S = 300000) {
  
  sigma1 <- sqrt(sigma1_sq)
  sigma2 <- sqrt(sigma2_sq)
  Sigma  <- matrix(c(sigma1^2, rho*sigma1*sigma2,
                    rho*sigma1*sigma2, sigma2^2),
                  nrow = 2, byrow = TRUE)
  
  # Theta = (Theta1, Theta2) ~ MVN((mu1,mu2), Sigma)
  theta_mat <- mvrnorm(n = S, mu = c(mu1,mu2), Sigma = Sigma)
  theta1 <- theta_mat[,1]
  theta2 <- theta_mat[,2]
  
  hurdle_probability <- sigma_fun(theta1)
  poisson_mean <- exp(theta2)
  
  #  Y1 = Z1 * (1 + N1)
  z1 <- rbinom(S, size = 1, prob = hurdle_probability)
  n1 <- rpois(S, lambda = poisson_mean)
  y1 <- z1 * (1 + n1)
  
  # Generate Y2 independently of Y1 conditional on the random effects
  z2 <- rbinom(S, size = 1, prob = hurdle_probability)
  n2 <- rpois(S, lambda = poisson_mean)
  y2 <- z2 * (1 + n2)
  
  return(conditional_summary(y1, y2))
}

rho_values <- c(-0.8, -0.5, 0, 0.5, 0.8)
results <- data.frame(
  rho = rho_values,
  mean_y2_given_y1_eq_1 = NA,
  mcse_y2_given_y1_eq_1 = NA,
  mean_y2_given_y1_eq_2 = NA,
  mcse_y2_given_y1_eq_2 = NA
)


set.seed(123)
for(i in seq_along(rho_values)) {
  out <- simulate_once(rho = rho_values[i])
  results$mean_y2_given_y1_eq_1[i] <- out["mean1"]
  results$mcse_y2_given_y1_eq_1[i] <- out["mcse1"]
  results$mean_y2_given_y1_eq_2[i] <- out["mean2"]
  results$mcse_y2_given_y1_eq_2[i] <- out["mcse2"]
}

print(round(results, 4))
