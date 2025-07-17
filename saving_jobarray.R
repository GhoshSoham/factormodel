source("Generate_subspace.R")
source("gibbs_subspace.R")
source("helpersampling.R")


# Posterior sampling and computing ratios
MCMCsub <- function(n, d, k, flag = 0, thin = 10, burnin = 1000, No_Iter = 10000) {
  # Generating samples
  dat <- GenerateData(n = n, d = d, k0 = k, r = 10, sigma = 0.01, flag = 0)

  intpara <- InitParam(d, k, lambda = 1)

  sampling <- factGibbsMod(
    X = dat$X, k = k, U0 = intpara$U, r0 = intpara$r,
    sigma0 = 0.01, Ut = dat$U, No_Iter = No_Iter,
    burnin = burnin, thin = thin
  )

  dist_summ <- c(
    mean(unlist(sampling$dist)), var(unlist(sampling$dist)), median(unlist(sampling$dist)),
    quantile(unlist(sampling$dist), names = FALSE, probs = c(0.05, 0.25, 0.75, 0.95))
  )

  r_summ <- c(
    mean(unlist(sampling$R)), var(unlist(sampling$R)), median(unlist(sampling$R)),
    quantile(unlist(sampling$R), names = FALSE, probs = c(0.05, 0.25, 0.75, 0.95))
  )


  return(rbind(dist_summ, r_summ, sampling$distpc, sampling$rest))
}


# your_script.R
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])

library(truncnorm)
library(vMF)
library(MASS)

# set.seed(task_id)
result <- MCMCsub(n = 100, d = 200, k = 1, flag = 0, thin = 10, burnin = 100, No_Iter = 1000)
colnames(result) <- c("mean", "var", "median", "0.05th", "0.25th", "0.75th", "0.95th")


# Append results to a central file (ensure it's thread-safe)
write.table(result,
  file = "n500d2000k1r10sig2.csv",
  append = TRUE,
  row.names = TRUE,
  col.names = !file.exists("n500d2000k1r10sig2.csv"), # Add headers only once
  sep = ","
)
