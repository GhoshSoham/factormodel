source('Generate_subspace.R')
source('gibbs_subspace.R')
source('helpersampling.R')

# parallel
library(foreach)
library(doParallel)

# Posterior sampling and computing ratios
MCMCsub <- function(n, d, k, flag = 0, thin = 10, burnin = 1000, No_Iter = 10000) {
  # Generating samples 
  dat = GenerateData(n = n, d = d, k0 = k, r = 5, sigma = 0.001, flag = 0)
  
  intpara = InitParam(d, k, lambda = 1)
  
  sampling = factGibbsMod(X = dat$X, k = k, U0 = intpara$U, r0 = intpara$r,
                          sigma0 = 0.01, Ut = dat$U, No_Iter = No_Iter, 
                          burnin = burnin, thin = thin)
  
  dist_summ = c(mean(unlist(sampling$dist)), var(unlist(sampling$dist)), median(unlist(sampling$dist)), 
              quantile(unlist(sampling$dist), names = FALSE, probs = c(0.05, 0.25, 0.75, 0.95)))
  
  r_summ = c(mean(unlist(sampling$R)), var(unlist(sampling$R)), median(unlist(sampling$R)), 
                quantile(unlist(sampling$R), names = FALSE, probs = c(0.05, 0.25, 0.75, 0.95)))
  
  distcom_summ = c(mean(unlist(sampling$distcom)), var(unlist(sampling$distcom)), median(unlist(sampling$distcom)), 
                quantile(unlist(sampling$distcom), names = FALSE, probs = c(0.05, 0.25, 0.75, 0.95)))
  
  return(rbind(dist_summ, r_summ, distcom_summ))
}

library(truncnorm)
library(vMF)
library(MASS)
result = MCMCsub(n = 50, d = 2000, k = 2, flag = 0, thin = 1, burnin = 0, No_Iter = 200)

# Execute the code

nworkers <- detectCores() - 1 # extra trick if running on your m 
cl <- makePSOCKcluster(nworkers)
registerDoParallel(cl)


nrep <- 50

output1 <- foreach (rep = 1:nrep) %dopar% {
  library(truncnorm)
  library(vMF)
  library(MASS)
  
  # Doing replicated simulation study
  MCMCsub(n = 100, d = 5000, k = 2, flag = 0, thin = 10, burnin = 1000, No_Iter = 10000)
}


re_dist = matrix(0, nrep, 7)
re_r = matrix(0, nrep, 7)
for(i in 1:nrep) {
  re_dist[i, ] = output1[[i]][1, ]
  re_r[i, ] = output1[[i]][2, ]
}
re_dist = round(re_dist, 6)
re_r = round(re_r, 6)

re_dist = as.data.frame(re_dist)
re_r = as.data.frame(re_r)

colnames(re_dist) = c("mean", "var", "median", "0.05th", "0.25th", "0.75th", "0.95th")
colnames(re_r) = c("mean", "var", "median", "0.05th", "0.25th", "0.75th", "0.95th")


write.csv(re_dist, file = "n100d5000k2dist.csv")
write.csv(re_r, file = "n100d5000k2r.csv")

nrep <- 50

output1 <- foreach (rep = 1:nrep) %dopar% {
  library(truncnorm)
  library(vMF)
  library(MASS)
  
  # Doing replicated simulation study
  MCMCsub(n = 200, d = 5000, k = 2, flag = 0, thin = 10, burnin = 1000, No_Iter = 10000)
}


re_dist = matrix(0, nrep, 7)
re_r = matrix(0, nrep, 7)
for(i in 1:nrep) {
  re_dist[i, ] = output1[[i]][1, ]
  re_r[i, ] = output1[[i]][2, ]
}
re_dist = round(re_dist, 6)
re_r = round(re_r, 6)

re_dist = as.data.frame(re_dist)
re_r = as.data.frame(re_r)

colnames(re_dist) = c("mean", "var", "median", "0.05th", "0.25th", "0.75th", "0.95th")
colnames(re_r) = c("mean", "var", "median", "0.05th", "0.25th", "0.75th", "0.95th")


write.csv(re_dist, file = "n200d5000k2dist.csv")
write.csv(re_r, file = "n200d5000k2r.csv")

nrep <- 50

output1 <- foreach (rep = 1:nrep) %dopar% {
  library(truncnorm)
  library(vMF)
  library(MASS)
  
  # Doing replicated simulation study
  MCMCsub(n = 500, d = 5000, k = 2, flag = 0, thin = 10, burnin = 1000, No_Iter = 10000)
}


re_dist = matrix(0, nrep, 7)
re_r = matrix(0, nrep, 7)
for(i in 1:nrep) {
  re_dist[i, ] = output1[[i]][1, ]
  re_r[i, ] = output1[[i]][2, ]
}
re_dist = round(re_dist, 6)
re_r = round(re_r, 6)

re_dist = as.data.frame(re_dist)
re_r = as.data.frame(re_r)

colnames(re_dist) = c("mean", "var", "median", "0.05th", "0.25th", "0.75th", "0.95th")
colnames(re_r) = c("mean", "var", "median", "0.05th", "0.25th", "0.75th", "0.95th")


write.csv(re_dist, file = "n500d5000k2dist.csv")
write.csv(re_r, file = "n500d5000k2r.csv")

stopCluster(cl) # Stop the cluster
