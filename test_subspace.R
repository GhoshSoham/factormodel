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
  distcom <- list()
  distlog <- list()
  
  U_cur <- U0
  r_cur <- r0
  sigma_cur = sigma0
  
  df_origin = scale(X, scale = FALSE)
  pc_sub = svd(df_origin, 0, k)$v
  
  for (iter in 1:(No_Iter + burnin + 1)) {
    # Update of latent variables G
    rsig = r_cur^2 + sigma_cur^2
    G_temp = sigma_cur/sqrt(rsig)*rnorm(n*k, 0, 1)
    G_temp = matrix(G_temp, nrow = n, ncol = k)
    datorth = X %*% U_cur
    G_cur = G_temp + datorth*(r_cur/rsig)  
    
    facttrace = sum(G_cur^2)
    factdata = sum(G_cur*datorth)
    
    
    # Update of r
    V2 = sigma_cur^2/facttrace
    mu2 = (factdata - sigma_cur^2)/facttrace
    r_cur <- rtruncnorm(1, mean = mu2, sd = V2,  a = sigma_cur)
    
    # Update of sigma^2
    # alpha = (d*n)/2 - 1 # shape param
    # beta = (sum(X^2) - 2*r_cur*factdata + r_cur^2*facttrace)/2 #scale param
    # a = 1/min(r_cur^2, 2/d)
    # b = d
    # sigma_temp = myrtgamma(n = 1, shape = alpha, scale = beta, lb = a, ub = b)
    # sigma_cur = 1/sigma_temp
    
    # Update of U using gibbs sampling
    M = (r_cur/sigma_cur^2)*t(X) %*% G_cur
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
      logpost <- - n*k*(log(rsig)/2) + (r_cur^2/rsig)*(sum(datorth^2)/(2*sigma_cur^2)) - r_cur
      loglik <- - n*k*(log(r0^2 + sigma0^2)/2) + 
        (r0^2/(r0^2 + sigma0^2))*(sum((X %*% U0)^2)/(2*sigma0^2))
      # Calculating distance between the projections
      # dist[[length(R) + 1]] <- (k - norm(crossprod(U_cur, Ut), type = "F")^2)
      dist[[length(R) + 1]] <- k - sum(crossprod(U_cur, Ut)^2)
      distcom[[length(R) + 1]] <- abs(sum(crossprod(pc_sub, Ut)^2) - 
                                        sum(crossprod(U_cur, Ut)^2))
      distlog[[length(R) + 1]] <- abs(logpost - loglik)
    }
  }
  
  output <- list(U = U, R = R, dist = unlist(dist), distcom = unlist(distcom), 
                 distlog = unlist(distlog) - mean(unlist(distlog)))
  
  return(output)
}

