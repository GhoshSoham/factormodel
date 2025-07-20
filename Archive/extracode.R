dense <- function(r, sigma, A) {
  rsig = r^2 + sigma^2 
  return(rsig^(- (k*n)/2)*exp(((r^2/rsig)*sum(A^2))/(2*sigma^2) - r))
}

## Gridded sampling
rsampling <- function(sigma, A) {
  rseq <- seq(1, 1000, length.out = 100)
  prob <- lapply(rseq, dense, sigma = 0.01, A = datorth)
}

dense(r = 10, sigma = 0.01, A = datorth)


U_cur <- U0
# Given constants
sigma2 <- 0.01^2    # Variance term (sigma^2)
k <- 2         # Constant multiplier in the exponent
n <- 500        # Number of observations
a <- X %*% U_cur         # Parameter for the exponential term


# Define the posterior function (unnormalized)
posterior <- function(r) {
  if (r^2 + sigma2 <= 0) return(0)  # Handle invalid cases
  return((r^2 + sigma2)^(-(k * n) / 2) * 
           exp(-(1 / (2 * (r^2 + sigma2))) * a))
}

# Metropolis-Hastings function
metropolis_hastings <- function(iter = 10000, init = 0, proposal_sd = 1) {
  samples <- numeric(iter)
  samples[1] <- init
  for (i in 2:iter) {
    # Propose a new sample
    proposal <- rnorm(1, mean = samples[i - 1], sd = proposal_sd)
    
    # Compute acceptance ratio
    alpha <- min(1, posterior(proposal) / posterior(samples[i - 1]))
    
    # Accept or reject
    if (runif(1) < alpha) {
      samples[i] <- proposal
    } else {
      samples[i] <- samples[i - 1]
    }
  }
  return(samples)
}

# Run the sampler
set.seed(123)
mh_samples <- metropolis_hastings(iter = 10000, init = 0, proposal_sd = 0.5)

# Discard burn-in and plot results
burn_in <- 1000
final_samples <- mh_samples[-(1:burn_in)]

# Plot histogram of samples
hist(final_samples, breaks = 50, probability = TRUE, col = "lightblue",
     main = "Samples from Posterior", xlab = "r")

# Compare with posterior curve (optional)
r_grid <- seq(min(final_samples), max(final_samples), length.out = 1000)
unnormalized_posterior <- sapply(r_grid, posterior)
lines(r_grid, unnormalized_posterior / sum(unnormalized_posterior), col = "blue", lwd = 2)



# Constants
sigma2 <- (0.01)^2    # Variance term (sigma^2)
k <- 2         # Constant multiplier in the exponent
n <- 500        # Number of observations
datorth = X %*% Ut
a <- sum(datorth^2)         # Parameter for the exponential term


# Define the grid for r
r_grid <- seq(0.2, 10, length.out = 1000)  # Grid points

# Compute the log-posterior values
log_posterior <- -((k * n) / 2) * log(r_grid^2 + sigma2) - (a / (2 * (r_grid^2 + sigma2))) - r_grid

# Find the maximum log-posterior value
max_log_posterior <- max(log_posterior)

# Scale the log-posterior values for stability
#log_posterior_scaled <- log_posterior - max_log_posterior
log_posterior_scaled <- log_posterior/mean(log_posterior)

# Convert back to unnormalized posterior
unnormalized_posterior <- exp(log_posterior_scaled)

# Normalize the posterior
normalized_posterior <- unnormalized_posterior / sum(unnormalized_posterior)

# Plot the normalized posterior
plot(r_grid, normalized_posterior, type = "l", col = "blue", lwd = 2,
     main = "Scaled Posterior Distribution", xlab = "r", ylab = "Density")



#########
## r sampling algorithm
rsampling <- function(n, k, sigma, a) {
  # Define the grid for r
  r_grid <- seq(0.25, 10, length.out = 1000)  # Grid points
  
  # Compute the log-posterior values
  log_posterior <- -((k * n) / 2) * log(r_grid^2 + sigma^2) - (a / (2 * (r_grid^2 + sigma^2))) - r_grid
  
  # Find the maximum log-posterior value
  #max_log_posterior <- max(log_posterior)
  
  # Scale the log-posterior values for stability
  log_posterior_scaled <- log_posterior/mean(log_posterior)
  
  # Convert back to unnormalized posterior
  unnormalized_posterior <- exp(log_posterior_scaled)
  
  # Normalize the posterior
  normalized_posterior <- unnormalized_posterior / sum(unnormalized_posterior)
  
  return(sample(r_grid, size = 1, prob = normalized_posterior, replace = TRUE))
}



### log sum exponnetial trick 

# Constants
sigma2 <- 1    # Variance term (sigma^2)
k <- 5         # Constant multiplier in the exponent
n <- 10        # Number of observations
a <- 2         # Parameter for the exponential term

# Define the grid for r
r_grid <- seq(-5, 5, length.out = 1000)  # Grid points

# Compute the log-posterior values
log_posterior <- -((k * n) / 2) * log(r_grid^2 + sigma2) - (a / (2 * (r_grid^2 + sigma2)))

# Apply the Log-Sum-Exp trick
max_log_posterior <- max(log_posterior)  # Maximum log value for stability
log_sum_exp <- max_log_posterior + log(sum(exp(log_posterior - max_log_posterior)))

# Normalize the posterior using the LSE trick
normalized_log_posterior <- log_posterior - log_sum_exp
normalized_posterior <- exp(normalized_log_posterior)

# Plot the normalized posterior
plot(r_grid, normalized_posterior, type = "l", col = "blue", lwd = 2,
     main = "Posterior Distribution (LSE Trick)", xlab = "r", ylab = "Density")


#####



