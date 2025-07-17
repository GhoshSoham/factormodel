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

# result = NULL
# for(i in 1:3) {
#   q = q0 + i
#   
#   # Initializing the parameters
#   beta_init = InitParam(d, q)$beta
#   Sigma_init = InitParam(d, q)$sigma
#   
#   ## this gives us a net: No.Iter/thin posterior samples
#   ## Gibbs sampler for the regular dataset: y
#   obj = factGibbsMod(y = y, k = q, beta_init, Sigma_init, nu, C0, s, No.Iter, burnin, thin)
#   
#   # Getting ratios corresponding to the posterior samples of loading matrix after thinning 
#   singratio = SingularRatio(obj, d, q0, q)
#   
#   summary = c(mean(singratio), var(singratio), median(singratio), quantile(singratio, names = FALSE, probs = c(0.05, 0.25, 0.75, 0.95)))
#   
#   result = c(result, summary)
# }
