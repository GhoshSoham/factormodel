##############
# Function that generates data for factor model
# n - No of samples in data
# d - Dimension of the data
# q0 - True no of factor
# flag - Deciding homoscedastic or heteroscedastic error 
GenerateData <- function(n, d, q0, flag = 1) {
  if (flag == 0) {
    psi = rep(10, d)
  }
  if (flag == 1) {
    psi = (10:(10+d-1)/4)
  }
  
  # Generating loading matrix
  W_true = matrix(rnorm(d*q0, 0, 1), nrow = d, ncol = q0)
  
  # Generating data matrix
  y = matrix(0, nrow = n, ncol = d)
  for(i in 1:n){
    f = rnorm(q0, mean = 0, sd = 1)
    y[i,] = rmvnorm(n = 1, mean = W_true %*% f, sigma = diag(psi))
  }
  
  return(y)
}

#############
# Function that initialize parameters for Gibbs sampling
# d - Dimension of the data
# q - Assumed no of factor
InitParam <- function(d, q) {
  # init value for beta0 used in Gibbs sampler
  beta_init = matrix(0, nr = d, nc = q) 
  for(i in 1:q) {
    if(i == 1){
      beta_init[1, 1] = rtnorm(n = 1, mean = 0, sd =1, lower= 0)
    } else {
      lower.vec = c(rep(-Inf, i - 1), 0)
      beta_init[i, 1:i] = rtmvnorm(n = 1, mean = rep(0, i), sigma = diag(i),
                                   lower = lower.vec)
    }
  }
  for(i in (q+1):d) {
    beta_init[i,] = rmvnorm(n = 1, mean = rep(0, q), sigma = diag(q))
  }
  
  Sigma_init = rep(1, d)
  
  return(list(beta = beta_init, sigma = Sigma_init))
}
