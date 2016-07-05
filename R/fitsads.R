library(sads)

#' @export

fitsads <- function(x){         #need to generalize, so any dist can be used, i.e.
  x <- x[!x==0]                 # function should accept a list of distribution to fit,
  geom <- fitsad(x, "geom")     # or as a default, 1 or 2 default distributions
  ls <- fitsad(x, "ls")
  poilog <- fitsad(x, "poilog")
  lnorm <- fitsad(x, "lnorm")
  res <- list(geom, ls, poilog, lnorm)
  
} 

#' @export

AIC.sads <- function(x){
  dat = x
  names.x = names(x)
  len = length(x)
  mat = matrix(ncol = 4, nrow = len)
  rownames(mat) <- names.x
  m1 <- dat[[1]][[1]]@sad
  m2 <- dat[[1]][[2]]@sad
  m3 <- dat[[1]][[3]]@sad
  m4 <- dat[[1]][[4]]@sad
  colnames(mat) <- c(m1,m2,m3,m4)
  for (i in 1:len){
    s1 <- AIC(x[[i]][[1]])
    s2 <- AIC(x[[i]][[2]])
    s3 <- AIC(x[[i]][[3]])
    s4 <- AIC(x[[i]][[4]])
    mat[i,1] <- s1
    mat[i,2] <- s2
    mat[i,3] <- s3
    mat[i,4] <- s4
  }
  return(mat)
}  

#' @export

min.AIC <- function(x){
  t(sapply(seq(nrow(x)), function(i) {
    j <- which.min(x[i,])
    c(paste(rownames(x)[i], colnames(x)[j], sep='_'), x[i,j])
  }))
}

#' @export

sads.coefs <- function(x){
  dat = x
  names.x = names(x)
  len = length(x)
  mat = matrix(ncol = 9, nrow = len)
  rownames(mat) <- names.x
  colnames(mat) <- c("Richness","N","geom_prob","ls_N","ls_alpha",
                     "poilog_mu","poilog_sig","lnorm_meanlog","lnorm_sdlog")
  for (i in 1:len){
    Richness <- length(x[[i]][[1]]@data$x)
    N = x[[i]][[2]]@fullcoef[1]
    m1.prob <- x[[i]][[1]]@fullcoef[1]
    m2.N <- x[[i]][[2]]@fullcoef[1]
    m2.alpha <- x[[i]][[2]]@fullcoef[2]
    m3.mu <- x[[i]][[3]]@fullcoef[1]
    m3.sig <- x[[i]][[3]]@fullcoef[2]
    m4.meanlog <- x[[i]][[4]]@fullcoef[1]
    m4.sdlog <- x[[i]][[4]]@fullcoef[2]
    mat[i,1] <- Richness
    mat[i,2] <- N
    mat[i,3] <- m1.prob
    mat[i,4] <- m2.N
    mat[i,5] <- m2.alpha
    mat[i,6] <- m3.mu
    mat[i,7] <- m3.sig
    mat[i,8] <- m4.meanlog
    mat[i,9] <- m4.sdlog
    
  }
  mat <- as.data.frame(mat)
  return(mat)
} 