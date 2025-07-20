library(foreach)
library(doParallel)


# Rue's sampler for MVN
# sample below is a sample from N_d(Q^{-1}b, Q^{-1}) 
mvnrue = function(Q,b){ 
  d = length(b)
  L = t(chol(Q))                      # Q = LL'; R returns right Cholesky factor by default
  yr = solve(t(L),rnorm(d))           # z~N(0,I_d), L'y=z
  vr = solve(L,b)
  thetar = solve(t(L),vr)             # L'theta = v, where Lv = b
  samp = yr + thetar
  return(samp)
}


##############
# Function that do Gibbs sampling for factor model
# y - Data matrix
# k - No of factor assumed
# beta0 - Initial factor loading matrix

factGibbsMod <- function(X, k, U0, R0, Ut, Rt, lambda, sigma, No_Iter, burnin, thin) {
  
  # Getting dimension information
  n <- nrow(X)
  d <- ncol(X)
  
  ######
  U <- list()
  R <- list()
  dist <- list()
  
  U_cur <- U0
  R_cur <- diag(R0)
  
  
  for (iter in 1:(No_Iter + burnin + 1)) {
    # Update of latent variables F
    # V2 <- solve(diag(k) + (t(R_cur) %*% t(U_cur) %*% U_cur %*% R_cur) * (1 / sigma^2))
    # F_temp <- MASS::mvrnorm(n, rep(0, k), V2)
    # F_temp <- F_temp + t(V2 %*% t(R_cur) %*% t(U_cur) %*% t(X)) * (1 / sigma^2)
    # F_cur <- F_temp
    V2 = diag(k) + (t(R_cur) %*% t(U_cur) %*% U_cur %*% R_cur) * (1 / sigma^2) 
    mu2 = (X %*% U_cur %*% R_cur)* (1 / sigma^2)
    F_temp = apply(mu2, MARGIN = 1, mvnrue, Q = V2)
    F_cur = t(F_temp)
    
    
    # Update of U
    # V1 <- solve(d * diag(k) + (t(R_cur) %*% t(F_cur) %*% F_cur %*% R_cur) * (1 / sigma^2))
    # U_temp <- MASS::mvrnorm(d, rep(0, k), V1)
    # U_temp <- U_temp + t(V1 %*% t(R_cur) %*% t(F_cur) %*% X) * (1 / sigma^2)
    # U_cur <- U_temp
    V1 = d * diag(k) + (t(R_cur) %*% t(F_cur) %*% F_cur %*% R_cur) * (1 / sigma^2)
    mu1 = (t(X) %*% F_cur %*% R_cur)* (1 / sigma^2)
    U_temp = apply(mu1, MARGIN = 1, mvnrue, Q = V1)
    U_cur = t(U_temp)
    
    
    # Update of singular values R
    V3 <- solve((t(F_cur) %*% F_cur) * (t(U_cur) %*% U_cur) * (1 / sigma^2))
    mu3 <- V3 %*% diag((t(F_cur) %*% X %*% U_cur) * (1 / sigma^2) - lambda * diag(k))
    R_temp <- TruncatedNormal::rtmvnorm(1, mu = mu3, sigma = V3, lb = sigma * rep(1, k))
    R_cur <- diag(R_temp)
    
    if (iter > burnin && (iter - burnin) %% thin == 0) {
      U[[length(U) + 1]] <- U_temp
      R[[length(R) + 1]] <- R_temp
      dist[[length(R) + 1]] <- norm(U_cur %*% R_cur %*% t(U_cur) - Ut %*% diag(Rt) %*% t(Ut), type = "2")
    }
  }

  output <- list(U = U, R = R, dist = dist)
  
  return(output)
}

d0 = 50
sigma = 1/sqrt(d0)
params = init_param(d=d0, k=2, lambda = 1)
U0 = params$U
R0 = params$R

dat = GenerateData(n=3000, d=d0, k0=2, sigma, flag = 0)
te = (t(dat$X) %*% dat$X)/n

norm(cov(dat$X) - dat$U %*% diag(dat$R^2) %*% t(dat$U) - sigma^2 * diag(d0), type = "2")
result = factGibbsMod(dat$X, 2, U0, R0, dat$U, dat$R, lambda = 1, sigma=sigma, No_Iter=30000, burnin=2000, thin=100)

sum(X$U %*% t(X$U) - t(result$U[[250]]) %*% result$U[[250]])^2 
fun = function(A) {return(max(svd(A)$d))}
list = lapply(result$U, FUN = fun)
plot(1:500, list)


