##############
# Function that do Gibbs sampling for factor model
# y - Data matrix
# k - No of factor assumed
# beta0 - Initial factor loading matrix

factGibbsMod <- function(y, k, beta0, Sigma0, nu, C0, s, No.Iter, burnin, thin) {
  
  # Getting dimension information
  T <- nrow(y)
  m <- ncol(y)
  
  
  # Take care of the case when k = 0
  if (k == 0) {
    output <- list()
    output$beta <- NA
    output$F <- NA
    Sigma.temp <- numeric() # a draw from the posterior
    for (j in 1:m) {
      Sigma.temp[j] <- rigamma(
        n = 1, alpha = (nu + T) / 2,
        beta = (sum(y[, j]^2) + nu * s^2) / 2
      )
    }
    output$Sigma[[1]] <- Sigma.temp # split out a posterior draw for Sigma
    return(output)
  }
  
  
  ######
  # list objects for the simulated parameters initialized
  
  beta <- vector(mode = "list", length = No.Iter / thin)
  Sigma <- vector(mode = "list", length = No.Iter / thin)
  ##############################
  
  
  # set the current values of the parameters
  beta.cur <- beta0
  Sigma.cur <- Sigma0
  
  
  for (iter in 1:(No.Iter + burnin)) {
    # update F
    F.temp <- rmvnorm(
      n = T,
      mean = rep(0, k),
      sigma = solve(diag(k) + t(beta.cur) %*% diag(1 / Sigma.cur) %*% beta.cur)
    )
    
    
    F.temp <- F.temp + t(solve(diag(k) + t(beta.cur) %*% diag(1 / Sigma.cur) %*% beta.cur) %*% t(beta.cur) %*% diag(1 / Sigma.cur) %*% t(y))
    
    F.cur <- F.temp
    
    # update Sigma
    Sigma.temp <- numeric(m)
    
    for (i in 1:m) {
      if (i <= k) {
        Sigma.temp[i] <- rigamma(
          n = 1, alpha = (nu + T) / 2,
          beta = (nu * s^2 + crossprod(y[, i] - F.cur[, 1:i] %*% as.matrix(beta.cur[i, 1:i]))) / 2
        )
      } else {
        Sigma.temp[i] <- rigamma(
          n = 1, alpha = (nu + T) / 2,
          beta = (nu * s^2 + crossprod(y[, i] - F.cur %*% as.matrix(beta.cur[i, ]))) / 2
        )
      }
    }
    
    # update the Sigma list when necessary
    if (iter > burnin && (iter - burnin) %% thin == 1 && thin != 1) {
      Sigma[[(iter - burnin) %/% thin + 1]] <- Sigma.temp
    } else if (thin == 1 && (iter - burnin) > 0) {
      Sigma[[iter - burnin]] <- Sigma.temp
    }
    
    Sigma.cur <- Sigma.temp
    
    
    # update beta
    
    beta.temp <- matrix(0, nr = m, nc = k)
    
    for (i in 1:m) {
      if (i <= k) {
        Ci <- solve(1 / C0 * diag(i) + 1 / Sigma.cur[i] * crossprod(F.cur[, 1:i]))
        mi <- Ci %*% (1 / Sigma.cur[i] * t(F.cur[, 1:i]) %*% y[, i])
        
        if (k != 1) {
          beta.temp[i, i] <- ars(1, ldkernel,
                                 ldkernelprima,
                                 x = c(0.1, 500, 1000), m = length(c(0.1, 500, 1000)), lb = TRUE, ub = FALSE,
                                 xlb = 0, k = k, i = i, mu = mi[i], c = Ci[i, i]
          )
        } else {
          # to tackle the instability of ars when k = 1
          beta.temp[i, i] <- rtnorm(n = 1, mean = mi[i], sd = sqrt(Ci[i, i]), lower = 0)
        }
        
        if (i != 1) {
          beta.temp[i, 1:(i - 1)] <-
            rmvnorm(1,
                    mean = mi[1:(i - 1)] + (1 / Ci[i, i]) * (beta.temp[i, i] - mi[i]) * Ci[1:(i - 1), i],
                    sigma = Ci[1:(i - 1), 1:(i - 1)] - 1 / Ci[i, i] * tcrossprod(Ci[1:(i - 1), i], Ci[i, 1:(i - 1)])
            )
        }
      } else {
        Ci <- solve(1 / C0 * diag(k) + 1 / Sigma.cur[i] * crossprod(F.cur))
        mi <- Ci %*% (1 / Sigma.cur[i] * t(F.cur) %*% y[, i])
        beta.temp[i, 1:k] <- rmvnorm(n = 1, mean = as.numeric(mi), sigma = Ci)
      }
    }
    
    beta.cur <- beta.temp
    
    
    if (iter > burnin && (iter - burnin) %% thin == 1 && thin != 1) {
      beta[[(iter - burnin) %/% thin + 1]] <- beta.temp
    } else if (thin == 1 && (iter - burnin) > 0) {
      beta[[iter - burnin]] <- beta.temp
    }
    
    
    # if(iter%%1000 == 0 )
    # cat("You are done with the ",iter, "th iterations when k is ", k,
    #       append= T,
    #     sep = "\n"
    #       )
  }
  
  output <- list()
  output$beta <- beta
  output$Sigma <- Sigma
  
  return(output)
}

#############
# Function that computes ratio corresponding to the singular values 
# d - Dimension of the data
# q0 - True no of factor
# q - Assumed no of factor
SingularRatio <- function(obj, d, q0, q) {
  net = length(obj$beta)
  singratio = rep(0, net)
  for (i in 1:net) {
    beta = matrix(unlist(obj$beta[i]), d, q)
    sing = svd(beta)$d
    singratio[i] = sum(sing[(q0+1):q])/sum(sing[1:q])
  }
  return(singratio)
}

#############
ldkernel <- function(x, k, i, mu, c){
  return((k - i)*log(x) - 0.5* (x - mu)^2/c)
}

#############
ldkernelprima <- function(x, k, i, mu, c){
  return((k - i)/x - (x- mu)/c )
}