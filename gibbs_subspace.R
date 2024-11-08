##############
# Function that do Gibbs sampling for factor model
# X - Data matrix of dimension n x d
# k - No of factor assumed
# U0 - Initial factor loading matrix
# r0 - Initial scaling matrix
# sigma0 - Initial residual st dev

factGibbsMod <- function(X, k, U0, r0, sigma0, Ut, No_Iter, burnin, thin) {
  
  # Getting dimension information
  n <- nrow(X)
  d <- ncol(X)
  
  ######
  U <- list()
  R <- list()
  dist <- list()
  
  U_cur <- U0
  r_cur <- r0
  sigma_cur = sigma0
  
  for (iter in 1:(No_Iter + burnin + 1)) {
    # Update of latent variables G
    V1 = ((sigma_cur^2)/(r_cur^2 + sigma_cur^2))*diag(k)
    G_temp <- mvrnorm(n, rep(0, k), V1)
    datorth = X %*% U_cur
    G_cur = G_temp + datorth*(r_cur/(r_cur^2 + sigma_cur^2))  
    
    facttrace = sum(G_cur^2)
    factdata = sum(G_cur*datorth)
    
    
    # Update of r
    V2 = sigma_cur^2/facttrace
    mu2 = (factdata - sigma_cur^2)/facttrace
    r_cur <- rtruncnorm(1, mean = mu2, sd = V2,  a = sigma_cur)
    
    # Update of sigma^2
    alpha = d*n - 1 # shape param
    beta = (sum(X^2) - 2*r_cur*factdata + r_cur^2*facttrace)/2 #scale param
    a = 1/min(r_cur^2, 10/d)
    b = d/1
    sigma_temp = myrtgamma(n = 1, shape = alpha, scale = beta, lb = a, ub = b)
    sigma_cur = 1/sigma_temp

    # Update of U using gibbs sampling
    M = (r_cur/sigma_cur^2)*t(X) %*% G_cur
    # Direction denotes the mean direction parameter
    M0 = matrix(rnorm(nrow(M)*ncol(M)), nrow(M), ncol(M))
    # If k = 1
    if(k == 1){
      U_cur = matrix(rvMF(1, theta = M), ncol = 1)
    } else {
      U_cur = rlangevin.gibbs(M, M0)
    }
    
    
    if (iter > burnin && (iter - burnin) %% thin == 0) {
      U[[length(U) + 1]] <- U_cur
      R[[length(R) + 1]] <- r_cur
      # Calculating distance between the projections
      dist[[length(R) + 1]] <- norm(t(U_cur) %*% Ut, type = "F")^2
      #dist[[length(R) + 1]] <- norm(U_cur %*% t(U_cur) - Ut %*% t(Ut), type = "2")
      print(iter)
    }
  }
  
  output <- list(U = U, R = R, dist = dist)
  
  return(output)
}

