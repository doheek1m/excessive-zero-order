## ---------------------------------------------------------------------------
## Supplementary Table S.2
## Bivariate-random effect analogue of Model 4 (Poisson-hurdle mixture with the
## softplus activation function of Remark 1) under the condition in (14):
##   Y_t = Z_t * (1 + N_t),
##     Z_t | Theta ~ Ber(logistic(Theta1)),  N_t | Theta ~ Pois(softplus(Theta2)),
##     (Theta1, Theta2) ~ MVN((0,0), Sigma).
## ---------------------------------------------------------------------------
 
library(MASS)

# Sigmoid function
sigma_fun <- function(x) 1 / (1 + exp(-x))

# Softplus function: lambda(x) = log(1 + exp(x))
softplus_fun <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

# Conditional means and Monte Carlo standard errors
conditional_summary <- function(y1, y2) {
  y2_given_y1_eq_0 <- y2[y1 == 0]
  y2_given_y1_eq_1 <- y2[y1 == 1]
  
  mean0 <- mean(y2_given_y1_eq_0)
  mean1 <- mean(y2_given_y1_eq_1)
  mcse0 <- sd(y2_given_y1_eq_0) / sqrt(length(y2_given_y1_eq_0))
  mcse1 <- sd(y2_given_y1_eq_1) / sqrt(length(y2_given_y1_eq_1))
  
  return(c(mean0 = mean0, mcse0 = mcse0,
           mean1 = mean1, mcse1 = mcse1))
}

simulate_once <- function(sigma1_sq, sigma2_sq = 1,
                          rho = 0.5, N = 2000000) {

  sigma1 <- sqrt(sigma1_sq)
  sigma2 <- sqrt(sigma2_sq)
  Sigma <- matrix(c(sigma1^2, rho*sigma1*sigma2,
                    rho*sigma1*sigma2, sigma2^2),
                  nrow = 2, byrow = TRUE)
  
  # Theta = (Theta1, Theta2) ~ MVN((0,0), Sigma)
  theta_mat <- mvrnorm(n = N, mu = c(0,0), Sigma = Sigma)
  theta1 <- theta_mat[,1]
  theta2 <- theta_mat[,2]
  
  # Bivariate-RE analogue of Model 4 with c1 = c2 = d1 = d2 = 0
  hurdle_probability_1 <- sigma_fun(theta1)
  poisson_mean_1 <- softplus_fun(theta2)
  hurdle_probability_2 <- sigma_fun(theta1)
  poisson_mean_2 <- softplus_fun(theta2)
  
  # Y1 = Z1 * (1 + N1)
  z1 <- rbinom(N, size = 1, prob = hurdle_probability_1)
  n1 <- rpois(N, lambda = poisson_mean_1)
  y1 <- z1 * (1 + n1)
  
  # Generate Y2 independently of Y1 conditional on the random effects
  z2 <- rbinom(N, size = 1, prob = hurdle_probability_2)
  n2 <- rpois(N, lambda = poisson_mean_2)
  y2 <- z2 * (1 + n2)
  
  return(conditional_summary(y1, y2))
}

sigma1_sq_values <- c(2.00, 1.00, 0.10, 0.01)
results <- data.frame(
  sigma1_sq = sigma1_sq_values,
  mean_y2_given_y1_eq_0 = NA,
  mcse_y2_given_y1_eq_0 = NA,
  mean_y2_given_y1_eq_1 = NA,
  mcse_y2_given_y1_eq_1 = NA
)

set.seed(123) 
for(i in seq_along(sigma1_sq_values)) {
  out <- simulate_once(sigma1_sq = sigma1_sq_values[i])
  results$mean_y2_given_y1_eq_0[i] <- out["mean0"]
  results$mcse_y2_given_y1_eq_0[i] <- out["mcse0"]
  results$mean_y2_given_y1_eq_1[i] <- out["mean1"]
  results$mcse_y2_given_y1_eq_1[i] <- out["mcse1"]
}

print(round(results, 4))



