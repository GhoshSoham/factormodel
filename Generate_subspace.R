##############
# Function that generates data for factor model
# n - No of samples in data
# d - Dimension of the data
# k0 - True no of factor
# flag - Deciding homoscedastic or heteroscedastic error 

GenerateData <- function(n, d, k0, r, sigma, flag = 0) {
  if (flag == 0) {
    psi = rep(sigma^2, d)
  }
  if (flag == 1) {
    psi = sample(round(seq(from = 0.05, to = 0.1, length.out = d), 3), size = d)
  }
  
  
  U = matrix(rnorm(d*k0, 0, 1), nrow = d, ncol = k0)
  W_true = r*svd(U)$u
  
  # Generating data matrix
  F = matrix(rnorm(n*k0, mean = 0, sd = 1), n, k0) 
  X = F %*% t(W_true) + matrix(rnorm(n*d, mean = 0, sd = 1), n, d)
  
  eff_rank = (k0*r^2 + d*sigma^2)/(r^2 + sigma^2)
  return(list(X = X, R = r, U = U, eff_rank = eff_rank))
}

#############
# Function that initialize parameters for Gibbs sampling
# d - Dimension of the data
# k - Assumed no of factor

InitParam <- function(d, k, lambda) {
  # Initialize U_init matrix with zeros
  #U_init <- matrix(rnorm(d*k), nrow = d, ncol = k)
  
  U = matrix(rnorm(d*k), nrow = d, ncol = k)
  U_init = svd(U)$u
  
  # Initialize R_init matrix with 1s 
  r_init <- rexp(1, lambda)
  
  return(list(U = U_init, r = r_init))
}

