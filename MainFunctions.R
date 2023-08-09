source("Generate.R")
source("Methods.R")

# parallel
library(foreach)
library(doParallel)

# Posterior sampling and computing ratios
MCMCRatio <- function(n, d, q0, flag = 1, thin = 20, burnin = 2000, No.Iter = 6000) {
  # Fixing some nuisance parameters
  nu = 2.2
  C0 = 1
  s = sqrt(.1/2.2)
  
  # Generating samples 
  y = GenerateData(n, d, q0, flag)
  
  result = NULL
  for(i in 1:3) {
    q = q0 + i
    
    # Initializing the parameters
    beta_init = InitParam(d, q)$beta
    Sigma_init = InitParam(d, q)$sigma
    
    ## this gives us a net: No.Iter/thin posterior samples
    ## Gibbs sampler for the regular dataset: y
    obj = factGibbsMod(y = y, k = q, beta_init, Sigma_init, nu, C0, s, No.Iter, burnin, thin)
    
    # Getting ratios corresponding to the posterior samples of loading matrix after thinning 
    singratio = SingularRatio(obj, d, q0, q)
    
    summary = c(mean(singratio), var(singratio), median(singratio), quantile(singratio, names = FALSE, probs = c(0.05, 0.25, 0.75, 0.95)))
    
    result = c(result, summary)
  }
  
  return(result)
}


# Execute the code

nworkers <- detectCores() # Initialize the cluster 
cl <- makeCluster(nworkers)
registerDoParallel(cl)


nrep <- 50

output <- foreach (rep = 1:nrep) %dopar% {
  library(mvtnorm)
  library(tmvtnorm)
  library(msm)
  library(pscl)
  library(ars)
  
  # Doing replicated simulation study
  MCMCRatio(n = 5000, d = 20, q0 = 5, flag = 1, thin = 20, burnin = 2000, No.Iter = 12000)
}


re = matrix(0, nrep, 21)
for(i in 1:nrep) {
  re[i, ] = output[[i]]
}
re1 = round(re[, 1:7], 6)
re2 = round(re[, 8:14], 6)
re3 = round(re[, 15:21], 6)
re1 = as.data.frame(re1)
re2 = as.data.frame(re2)
re3 = as.data.frame(re3)
colnames(re1) = c("mean", "var", "median", "0.05th", "0.25th", "0.75th", "0.95th")
colnames(re2) = c("mean", "var", "median", "0.05th", "0.25th", "0.75th", "0.95th")
colnames(re3) = c("mean", "var", "median", "0.05th", "0.25th", "0.75th", "0.95th")


write.csv(re1, file = "n5000d20qt5q6.csv")
write.csv(re2, file = "n5000d20qt5q7.csv")
write.csv(re3, file = "n5000d20qt5q8.csv")

stopCluster(cl) # Stop the cluster
