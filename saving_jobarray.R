source('Generate_subspace.R')
source('gibbs_subspace.R')
source('helpersampling.R')


# Posterior sampling and computing ratios
MCMCsub <- function(n, d, k, flag = 0, thin = 10, burnin = 1000, No_Iter = 10000) {
  # Generating samples 
  dat = GenerateData(n = n, d = d, k0 = k, r = 1, sigma = 0.001, flag = 0)
  
  intpara = InitParam(d, k, lambda = 1)
  
  sampling = factGibbsMod(X = dat$X, k = k, U0 = intpara$U, r0 = intpara$r,
                          sigma0 = 0.001, Ut = dat$U, No_Iter = No_Iter, 
                          burnin = burnin, thin = thin)
  
  dist_summary = c(mean(unlist(sampling$dist)), var(unlist(sampling$dist)), median(unlist(sampling$dist)), 
                   quantile(unlist(sampling$dist), names = FALSE, probs = c(0.05, 0.25, 0.75, 0.95)))
  
  r_summary = c(mean(unlist(sampling$R)), var(unlist(sampling$R)), median(unlist(sampling$R)), 
                quantile(unlist(sampling$R), names = FALSE, probs = c(0.05, 0.25, 0.75, 0.95)))
  
  return(rbind(dist_summary, r_summary))
}


# # your_script.R
# args <- commandArgs(trailingOnly = TRUE)
# task_id <- as.integer(args[1])

library(truncnorm)
library(vMF)
library(MASS)


result <- MCMCsub(n = 100, d = 1000, k = 2, flag = 0, thin = 1, burnin = 0, No_Iter = 100)
colnames(result) = c("mean", "var", "median", "0.05th", "0.25th", "0.75th", "0.95th")

microbenchmark::microbenchmark(MCMCsub(n = 2000, d = 3000, k = 2, flag = 0, thin = 1, burnin = 0, No_Iter = 20),
                               times = 2)

# Append results to a central file (ensure it's thread-safe)
write.table(result, 
            file = "n100d5k2.csv", 
            append = TRUE, 
            row.names = FALSE, 
            col.names = !file.exists("n100d5k2.csv"), # Add headers only once
            sep = ",")
