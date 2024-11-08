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

## truncated gamma sampling
myrtgamma = function(n, shape, scale, lb, ub) {
  num_intervals <- 100 # You can increase this for finer intervals
  # Generate interval edges
  interval_edges <- seq(lb, ub, length.out = num_intervals + 1)
  # Calculate probabilities for each interval
  interval_probs <- sapply(1:num_intervals, function(i) {
    # Integrate the density over the interval
    # lower <- interval_edges[i]
    # upper <- interval_edges[i + 1]
    # integrate(continuous_density, lower, upper)$value
    pgamma(interval_edges[i + 1], shape = shape, scale = scale) - 
      pgamma(interval_edges[i], shape = shape, scale = scale)
  })
  # Normalize to create a probability distribution over intervals
  interval_probs <- interval_probs / sum(interval_probs)
  sampled_intervals <- sample(1:num_intervals, size = n, replace = TRUE, prob = interval_probs)
  
  # Generate final samples by uniformly sampling within chosen intervals
  samples <- sapply(sampled_intervals, function(i) {
    runif(1, interval_edges[i], interval_edges[i + 1])
  })
  
}
