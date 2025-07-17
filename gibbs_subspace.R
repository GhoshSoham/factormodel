##############
# Function that do Gibbs sampling for factor model
# X - Data matrix of dimension n x d
# k - No of factor assumed
# U0 - Initial factor loading matrix
# r0 - Initial scaling matrix
# sigma0 - Initial residual st dev
# lambda - Parameter for exponential

factGibbsMod <- function(X, k, U0, r0, sigma0, Ut, No_Iter, burnin, thin) {
  
  # Getting dimension information
  n <- nrow(X)
  d <- ncol(X)
  
  ######
  U <- list()
  R <- list()
  dist <- list()
  distcom <- list()
  distr <- list()
  #distlog <- list()
  
  U_cur <- U0
  r_cur <- r0
  sigma_cur = sigma0
  
  df_origin = scale(X, scale = FALSE)
  #samsig = crossprod(df_origin)
  pc_sub = svd(df_origin, 0, k)$v
  pc_r = svd(df_origin, 0, k)$d[1]/sqrt(n)
  
  for (iter in 1:(No_Iter + burnin + 1)) {
    # Update of latent variables G
    rsig = r_cur^2 + sigma_cur^2
    G_temp = sigma_cur/sqrt(rsig)*rnorm(n*k, 0, 1)
    G_temp = matrix(G_temp, nrow = n, ncol = k)
    datorth = X %*% U_cur
    G_cur = G_temp + datorth*(r_cur/rsig)  
    
    facttrace = sum(G_cur^2)
    factdata = sum(G_cur*datorth)
    
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
    
    # Update of r
    # V2 = sigma_cur^2/facttrace
    # mu2 = (factdata - sigma_cur^2)/facttrace
    # r_cur <- rtruncnorm(1, mean = mu2, sd = V2,  a = sigma_cur)
    
    r_cur <- rsampling(n, k, sigma_cur, sum(datorth^2))
    #r_est <- sqrt(sum(datorth^2)/(n*k) - sigma_cur^2) # MoM wo origin shift
    #r_est = sqrt((sum(df_origin^2)/n - d*sigma_cur^2)/k) # MoM another way
    #r_est <- sqrt(sum((df_origin %*% U_cur)^2)/(n*k) - sigma_cur^2) # MoM w origin shift
    
    
    if (iter > burnin && (iter - burnin) %% thin == 0) {
      U[[length(U) + 1]] <- U_cur
      R[[length(R) + 1]] <- r_cur
      
      dist[[length(R) + 1]] <- sum((tcrossprod(U_cur, U_cur) - tcrossprod(Ut, Ut))^2)
      
    }
  }
  distpc <- sum((tcrossprod(pc_sub, pc_sub) - tcrossprod(Ut, Ut))^2)
  output <- list(U = U, R = unlist(R), dist = unlist(dist), distpc = distpc,
                 rest = pc_r)
  
  return(output)
}

