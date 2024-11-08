library(truncnorm)
library(vMF)
library(MASS)

source('Generate_subspace.R')
source('gibbs_subspace.R')
source('helpersampling.R')

n = 100
d = 5000
k = 2

dat = GenerateData(n = n, d = d, k0 = 2, r = 5, sigma = 0.001, flag = 0)

intpara = InitParam(d, k, lambda = 1)

sampling = factGibbsMod(X = dat$X, k = k, U0 = intpara$U, r0 = intpara$r,
             sigma0 = 0.001, Ut = dat$U, No_Iter = 10, burnin = 5, thin = 2)

save(sampling, file = "/Users/sohamghosh/Desktop/factormodel/sampling.RData")
load("/Users/sohamghosh/Desktop/Cluster/sampling.RData")
