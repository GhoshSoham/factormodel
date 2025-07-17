## Von Mises sampling
rlangevin.gibbs <- function(M, X) {
  sM <- svd(M)
  H <- sM$u %*% diag(sM$d)
  Y <- X %*% sM$v
  
  m <- nrow(H)
  R <- ncol(H)
  for (r in sample(1:R)) {
    N <- MASS::Null(Y[, -r])
    y <- matrix(vMF::rvMF(size = 1, theta = t(N) %*% H[, r]), ncol = 1)
    Y[, r] <- N %*% y
  }
  return(Y %*% t(sM$v))
}

## r sampling algorithm
rsampling <- function(n, k, sigma, a) {
  # Define the rate parameter lambda
  # lambda <- 2  # Example rate parameter
  # nu = 3 # df for t distribution
  
  # Number of grid points
  N <- 2000
  
  # Generate uniform quantiles (exclude 1 to avoid infinity at the tail)
  uniform_quantiles <- seq(0, 1, length.out = N + 1)[-c(1, (N+1))]  # Avoid q = 1
  
  # Transform using the inverse CDF of the exponential distribution
  # Instead of uniform grid, fix the grid using quantiles of exponential
  # r_grid <- -(1 / lambda) * log(1 - uniform_quantiles)
  # Transform using the inverse CDF of the half-t distribution
  # r_grid <- qt(uniform_quantiles, df = nu)  # Quantile function of t-distribution
  # Transform using the inverse CDF of the half-Cauchy distribution
  # r_grid <- tan(pi * uniform_quantiles / 2)
  
  
  # Compute the CDF at the truncation point sigma
  F_sigma <- (1/pi) * atan(sigma) + 0.5
  
  # Transform using the inverse CDF of the truncated Cauchy
  r_grid <- tan(pi * (uniform_quantiles * (1 - F_sigma) + F_sigma - 0.5))
  
  # Ensure values are at least sigma
  r_grid <- r_grid[r_grid >= sigma]
  
  # Compute the log-posterior values
  # exponential prior
  # log_posterior <- -((k * n) / 2) * log(r_grid^2 + sigma^2) - 
  #   (a / (2 * (r_grid^2 + sigma^2))) - lambda * r_grid
  # half t prior
  # log_posterior <- -((k * n) / 2) * log(r_grid^2 + sigma^2) - 
  #      (a / (2 * (r_grid^2 + sigma^2))) - (nu + 1)/2 * log(1 + r_grid^2/nu)
  # cauchy prior
  # log_posterior <- -((k * n) / 2) * log(r_grid^2 + sigma^2) -
  #      (a / (2 * (r_grid^2 + sigma^2))) - log(1 + r_grid^2)
  log_posterior <- -((k * n) / 2) * log(r_grid^2 + sigma^2) -
          (a / (2 * (r_grid^2 + sigma^2)))
  
  # Improper prior
  # log_posterior <- -((k * n) / 2) * log(r_grid^2 + sigma^2) -
  #         (a / (2 * (r_grid^2 + sigma^2))) 
  
   
  # Apply the Log-Sum-Exp trick
  max_log_posterior <- max(log_posterior)  # Maximum log value for stability
  log_sum_exp <- max_log_posterior + log(sum(exp(log_posterior - max_log_posterior)))
  
  # Normalize the posterior using the LSE trick
  normalized_log_posterior <- log_posterior - log_sum_exp
  normalized_posterior <- exp(normalized_log_posterior)
  
  return(sample(r_grid, size = 1, prob = normalized_posterior, replace = TRUE))
}
