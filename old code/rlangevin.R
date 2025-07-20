#####
library(MASS)
#####


rW <- function(kap, m) {
  # Set RNG state for reproducibility
  # set.seed(123)  # You can choose your seed value

  # Calculate b, x0, and c
  b <- (-2.0 * kap + sqrt(4 * kap^2 + (m - 1)^2)) / (m - 1)
  x0 <- (1 - b) / (1 + b)
  c <- kap * x0 + (m - 1) * log(1 - x0^2)

  done <- FALSE
  W <- numeric(1) # Initialize W

  # Sampling loop
  while (!done) {
    Z <- rbeta(1, (m - 1) / 2, (m - 1) / 2) # Sample from beta distribution
    W[1] <- (1 - (1 + b) * Z) / (1 - (1 - b) * Z) # Calculate W[0]
    U <- runif(1) # Sample from uniform distribution

    if (kap * W[1] + (m - 1) * log(1 - x0 * W[1]) - c > log(U)) {
      done <- TRUE
    }
  }

  return(W)
}


rvmf <- function(kmu) {
  kap <- sqrt(sum(kmu^2))
  mu <- kmu / kap
  m <- length(mu)
  if (kap == 0) {
    u <- rnorm(length(kmu))
    u <- u / sqrt(sum(u^2))
  }
  if (kap > 0) {
    if (m == 1) {
      u <- (-1)^rbinom(1, 1, 1 / (1 + exp(2 * kap * mu)))
    }

    if (m > 1) {
      W <- rW(kap, m)
      V <- rnorm(m - 1)
      V <- V / sqrt(sum(V^2))
      x <- c((1 - W^2)^.5 * t(V), W)
      u <- cbind(MASS::Null(mu), mu) %*% x
    }
  }
  u
}
#####

#####
rlangevin <- function(M) {
  svdM <- svd(M)
  H <- svdM$u %*% diag(svdM$d)
  m <- dim(H)[1]
  R <- dim(H)[2]

  ##
  cmet <- FALSE
  rej <- 0
  while (!cmet) {
    U <- matrix(0, m, R)
    U[, 1] <- rvmf(H[, 1])


    lr <- 0

    for (j in 2:R) {
      N <- Null(U[, 1:(j - 1)])
      x <- rvmf(t(N) %*% H[, j])
      U[, j] <- N %*% x

      xn <- sqrt(sum((t(N) %*% H[, j])^2))
      xd <- sqrt(sum(H[, j]^2))
      lbr <- log(besselI(xn, .5 * (m - j - 1), expon.scaled = T)) -
        log(besselI(xd, .5 * (m - j - 1), expon.scaled = T))
      if (is.na(lbr)) {
        lbr <- .5 * (log(xd) - log(xn))
      }
      lr <- lr + lbr + (xn - xd) + .5 * (m - j - 1) * (log(xd) - log(xn))
    }

    cmet <- (log(runif(1)) < lr)
    rej <- rej + (1 - 1 * cmet)
  }
  list(U = U %*% t(svd(M)$v), rej = rej)
}
#####

C <- matrix(rnorm(40), nrow = 20, ncol = 2)
X <- matrix(rnorm(40), nrow = 20, ncol = 2)



rlangevin.gibbs.old <- function(M, X) {
  sM <- svd(M)
  H <- sM$u %*% diag(sM$d)
  Y <- X %*% sM$v

  m <- dim(H)[1]
  R <- dim(H)[2]
  for (r in sample(1:R)) {
    N <- MASS::Null(Y[, -r])
    y <- rvmf(t(N) %*% H[, r])
    Y[, r] <- N %*% y
  }
  Y %*% t(sM$v)
}

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

## Inv gamma sampling
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

samples = myrtgamma(n = 5000, shape = 2, scale = 3, lb = 1, ub = 10)
hist(samples, breaks = 50, probability = TRUE, main = "Sampled from Discretized Truncated Normal")
curve(dgamma(x, shape = 2, scale = 3), from = 1, to = 10, add = TRUE, col = "blue", lwd = 2)

