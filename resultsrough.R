library(truncnorm)
library(vMF)
library(MASS)

source('Generate_subspace.R')
source('gibbs_subspace.R')
source('helpersampling.R')

n = 100
d = 500
k = 1
lambda = 1

dat = GenerateData(n = n, d = d, k0 = 1, r = 1, sigma = 0.1, flag = 0)

intpara = InitParam(d, k, lambda)

sampling = factGibbsMod(X = dat$X, k = k, U0 = intpara$U, r0 = intpara$r,
                        sigma0 = 0.1, Ut = dat$U, No_Iter = 2000, burnin = 200, thin = 1)

## r = 1.5 r = 2, r = 2.5, r = 3

microbenchmark::microbenchmark(factGibbsMod(X = dat$X, k = k, U0 = intpara$U, r0 = intpara$r,
                                            sigma0 = 0.01, Ut = dat$U, No_Iter = 200, burnin = 0, thin = 1), times = 2)

save(sampling, file = "/Users/sohamghosh/Desktop/factormodel/sampling.RData")
get(load("/Users/sohamghosh/Desktop/Cluster/samplingdat.RData"))

X = dat$X; k = k; U0 = dat$U; r0 = dat$R;
sigma0 = 0.01; Ut = dat$U; No_Iter = 30; burnin = 10; thin = 1


te = prcomp(dat$X)

(norm(t(te3) %*% dat$U, type = "F")^2)

dat1 = scale(dat$X, scale = FALSE)
te3 = svd(dat1, 0, 2)$v

