##############
# Function that generates data for factor model
# n - No of samples in data
# d - Dimension of the data
# q0 - True no of factor
# flag - Deciding homoscedastic or heteroscedastic error 
GenerateData <- function(n, d, k0, sigma, flag = 0) {
  if (flag == 0) {
    psi = rep(sigma^2, d)
  }
  if (flag == 1) {
    psi = sample(round(seq(from = 0.05, to = 0.1, length.out = d), 3), size = d)
  }
  
  # U <- matrix(0, nrow = d, ncol = k0)
  # for (i in 1:d) {
  #   U[i, ] <- MASS::mvrnorm(1, rep(0, k0), diag(k0))
  # }
  # U = U/sqrt(d)
  # R = rexp(k0, 1)
  # W_true = U %*% diag(R)
  
  U = pracma::orth(matrix(rnorm(d*k0), d, k0))
  #R = rexp(k0, lambda)
  R = sample(10:20, k0, replace = TRUE)
  W_true = U %*% diag(R)
  
  # Generating data matrix
  F = matrix(rnorm(n*k0, mean = 0, sd = 1), n, k0) 
  X = F %*% t(W_true) + matrix(rnorm(n*d, mean = 0, sd = 1), n, d)
  
  eff_rank = (sum(R^2) + d*sigma^2)/(max(R^2) + sigma^2)
  return(list(X = X, R = R, U = U, eff_rank = eff_rank))
}

#############
# Function that initialize parameters for Gibbs sampling
# d - Dimension of the data
# q - Assumed no of factor
init_param <- function(d, k, lambda) {
  # Initialize U_init matrix with zeros
  U_init <- matrix(0, nrow = d, ncol = k)
  
  for (i in 1:d) {
    U_init[i, ] <- MASS::mvrnorm(1, rep(0, k), diag(k))
  }
  U_init = U_init/sqrt(d)
  
  # Initialize R_init matrix with 1s 
  R_init <- rexp(k, lambda)
  
  return(list(U = U_init, R = R_init))
}

