## Von Mises sampling
rlangevin.gibbs <- function(M, X) {
  sM <- svd(M)
  H <- sM$u %*% diag(sM$d)
  Y <- X %*% sM$v
  
  m <- nrow(H)
  R <- ncol(H)
  for (r in sample(1:R)) {
    N <- MASS::Null(Y[, -r])
    y <- matrix(vMF::rvMF(size = 1, theta = t(N) %*% H[, r]), ncol = 1)
    Y[, r] <- N %*% y
  }
  return(Y %*% t(sM$v))
}

