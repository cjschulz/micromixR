fit_tnorm2 <- function(x) {
  x=x
  xmin = if(min(x) < log2(0.5)){
    min(x) -0.5
  } else {
    xmin = log2(0.5)
  }
  dtnorm0 <- function(X, mean, sd, log = FALSE) {dtnorm(X, mean, sd, xmin, Inf, log)} 
  res1 <- tryCatch(fitdistr(x, dtnorm0, start=list(mean=mean(x), sd=sd(x))),
                          error=function(e) fitdist(x, dnorm))
  res1$estimate[3] <- xmin
  res1$estimate[4] <- Inf
  names(res1$estimate)[3:4] <- c("lower", "upper")
  res1$sd[3] <- "NA"
  res1$sd[4] <- "NA"
  return(res1)
} #needs to add lower and upper params to list

fit.tNorm <- function(x) apply(x, 2, fit_tnorm2)



#rmvnorm from Spieceasi package
rmvnorm.se <- function(n=100, mu=rep(0,10), Sigma=diag(10), tol=1e-6, empirical=TRUE) {
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p))) 
    stop("incompatible arguments")
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  X <- matrix(rnorm(p * n), n)
  if (empirical) {
    X <- scale(X, TRUE, FALSE)
    X <- X %*% svd(X, nu = 0, nv = length(mu))$v
    X <- scale(X, FALSE, TRUE)
  }
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
  return(t(X))
}



