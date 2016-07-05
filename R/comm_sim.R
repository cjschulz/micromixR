### microbial community simulator

trim_otus <- function (x){
  x[ rowSums(x)!=0, colSums(x)!=0 ]
  x[ ,colSums(x==0) <= nrow(x)-2]
}


###need method to handle zeros in colSums###

fit_comm <- function(x, normalize=T) {
  mat <- x
  x.small <- min(rowSums(mat))
  if(normalize){
  x.norm <- as.data.frame(t(apply(mat, 1, function(x) rmultinom(1, size = x.small, prob = x))))
  colnames(x.norm) <- colnames(mat)
  } else { x.norm <- mat
  }
  x.norm <- trim_otus(x.norm)
  norm_otus <- list(x.norm)
  class(norm_otus) <- "otu_table"
  fitNB_x <- fitNB.mle(x.norm)
  fitPois_x <- fitPois.mle(x.norm)
  fitLN_x <- fitLN.mle(x.norm)
  fitDM_x <- fitDM(x.norm)
  fitCLR_x <- fitCLR(mat) #need to add method
  res <- list(norm_otus, fitNB_x, fitPois_x, fitLN_x, fitDM_x, fitCLR_x)
  return(res)
}

AIC.fit_comm <- function(x){
  dat = x
  names.x = names(x[[2]])
  len = length(x[[2]])
  mat = matrix(ncol = 4, nrow = len)
  rownames(mat) <- names.x
  m1 <- dat[[2]][[1]]$distname
  m2 <- dat[[3]][[2]]$distname
  m3 <- dat[[4]][[3]]$distname
  m4 <- dat[[6]][[4]]$distname
  colnames(mat) <- c(m1,m2,m3,m4)
  for (i in 1:len){
    s1 <- x[[2]][[i]]$aic
    s2 <- x[[3]][[i]]$aic
    s3 <- x[[4]][[i]]$aic
    s4 <- x[[6]][[i]]$aic
    mat[i,1] <- s1
    mat[i,2] <- s2
    mat[i,3] <- s3
    mat[i,4] <- s4
  }

  return(mat)
}

# fit_comm1 <- function(x) { #works, but not normalized
#   fitNB_x <- fitNB.mle(x)
#   fitPois_x <- fitPois.mle(x)
#   fitLN_x <- fitLN.mle(x)
#   fitDM_x <- fitDM(x)
#   res <- list(fitNB_x, fitPois_x, fitLN_x, fitDM_x)
#   return(res)
# }


fitNB.mme <- function(x) {
  x[ rowSums(x)!=0, colSums(x)!=0 ]
  fitneg <- function(x) {
  require(fitdistrplus)
  fitdist(x, "nbinom", method="mme")
  }
  fit <- apply(x, 2, fitneg)
  class(fit) <- c("comm.params","nbinom")
  return(fit)
}

fitNB.mle <- function(x) {
  x[ rowSums(x)!=0, colSums(x)!=0 ]
  fitneg <- function(x) {
  require(fitdistrplus)
  fitdist(x, "nbinom", method="mle")
  }
  fit <- apply(x, 2, fitneg)
  class(fit) <- c("comm.params","nbinom")
  return(fit)
}


fitPois.mle <- function(x) {
  x[ rowSums(x)!=0, colSums(x)!=0 ]
  fitp <- function(x) {
  require(fitdistrplus)
  fitdist(x, "pois", method="mle")
  }
  fit <- apply(x, 2, fitp)
  class(fit) <- c("comm.params","Pois")
  return(fit)
}


# fit normal distribution, for use with clr transformed data
fitCLR <- function(x) {
  xCLR <- clr_samples_rows(x)
  fitC <- fitNorm.mle(xCLR)
  class(fitC) <- c("comm.params","CLR")
  return(fitC)
}



fitNorm.mle <- function(x) {
  fitN <- function(x) {
  require(fitdistrplus)
  fitdist(x, "norm", method="mle")
  }
  fit <- apply(x, 2, fitN)
  class(fit) <- c("comm.params","Norm")
  return(fit)
}


fitLN.mle <- function(x) {
  x <- x[ rowSums(x)!=0, colSums(x)!=0 ]
  x <- x+0.5
  fitl <- function(x) {
    require(fitdistrplus)
    fitdist(x, "lnorm", method="mle")
  }
  fit <- apply(x, 2, fitl)
  class(fit) <- c("comm.params","LN")
  return(fit)
}

fitDM <- function(x) {
  x[ rowSums(x)!=0, colSums(x)!=0 ]
  require(HMP)
  res <- DM.MoM(x)
  class(res) <- c("comm.params", "DM")
  return(res)
}

#Fits a log-normal distribution
#fitLN <- function(x) {
#  x[ rowSums(x)!=0, colSums(x)!=0 ]
#  fitL <- function(x) {
#    require(MASS)
#    fitdistr(x, "log-normal")
#  }
#  fit <- apply(x, 2, fitL)
#}


### simulate a microbial community

# takes results from fitted distribution(s), for example a list from fitNB
# returns a simulated community
# # samples can be changed
# each simulated community can be named for tracking later

# simulate_comm_nbinom <- function(comm_params, samples=100, ids="sim_comm") {
#   paramat <- do.call('rbind', comm_params)
#   paramat2 <- as.list(paramat[,1])
#   paramat3 <- as.data.frame(do.call('rbind', paramat2))
#   size <- as.list(paramat3$size)
#   mu <- as.list(paramat3$mu)
#   test1 <- lapply(seq_along(mu), function(i) {
#     rnbinom(samples, mu=mu[[i]], size=size[[i]])
#   })
#   test1.results <- do.call('rbind', test1)
#   test1.results <- as.data.frame(t(test1.results))
#   names(test1.results) <- rownames(paramat3)
#   newnames <- paste(ids, seq(1:nrow(test1.results)), sep="")
#   row.names(test1.results) <- newnames
#   results <- test1.results
#
# }


### create synthetic microbial samples based given the following:

# richness (# otus)
# number of samples
# depth (otu counts per sample)
# either exponential or power law parameters for rank abundance simulation
# size parameters (i.e., lognormal simulation)



## Simulate communities using distributions of paramters for NB

simulate_comm_sd <- function(comm_params, samples=100, ids="sim_comm") {
  require(truncnorm)
  test1 <- lapply(seq_along(comm_params), function(i) {
    rnbinom(samples, mu=comm_params[[i]][[1]][[2]],
            size=(rtruncnorm(1, a=0, b=Inf, mean=comm_params[[i]][[1]][[1]], sd=comm_params[[i]][[2]][[1]])))
  })
  test1.results <- do.call('rbind', test1)
  test1.results <- as.data.frame(t(test1.results))
  names(test1.results) <- names(comm_params)
  newnames <- paste(ids, seq(1:nrow(test1.results)), sep="")
  row.names(test1.results) <- newnames
  results <- test1.results
}

# Methods for simulating communities using fitted paramters, or any list of suitable paramaters

simulate_comm <- function (x, ...) {
  UseMethod("simulate_comm", x)
}


simulate_comm.otu_table <- function(x, ...) {
  return(x[[1]])
}

simulate_comm.DM <- function(comm_params, samples=100, depth=1000, ids="DM_sim_comm", ...) {
  require(HMP)
  Nrs <- rep(depth, samples)
  test1 <- Dirichlet.multinomial(Nrs, comm_params[[2]])
  test1.results <- as.data.frame(test1)
  names(test1.results) <- names(comm_params[[2]])
  newnames <- paste(ids, seq(1:nrow(test1.results)), sep="")
  row.names(test1.results) <- newnames
  results <- test1.results
}

simulate_comm.nbinom <- function(comm_params, samples=100, ids="NB_sim_comm", ...) {
  test1 <- lapply(seq_along(comm_params), function(i) {
    rnbinom(samples, mu=comm_params[[i]][[1]][[2]],
            size=comm_params[[i]][[1]][[1]])
  })
  test1.results <- do.call('rbind', test1)
  test1.results <- as.data.frame(t(test1.results))
  names(test1.results) <- names(comm_params)
  newnames <- paste(ids, seq(1:nrow(test1.results)), sep="")
  row.names(test1.results) <- newnames
  results <- test1.results
}

simulate_comm.Pois <- function(comm_params, samples=100, ids="Pois_sim_comm", ...) {
  test1 <- lapply(seq_along(comm_params), function(i) {
    rpois(samples, lambda=comm_params[[i]][[1]][[1]])
  })
  test1.results <- do.call('rbind', test1)
  test1.results <- as.data.frame(t(test1.results))
  names(test1.results) <- names(comm_params)
  newnames <- paste(ids, seq(1:nrow(test1.results)), sep="")
  row.names(test1.results) <- newnames
  results <- test1.results
}

simulate_comm.LN <- function(comm_params, samples=100, ids="LN_sim_comm", ...) {
  test1 <- lapply(seq_along(comm_params), function(i) {
    rlnorm(samples, meanlog=comm_params[[i]][[1]][[1]],
           sdlog=comm_params[[i]][[1]][[2]])
  })
  test1.results <- do.call('rbind', test1)
  test1.results <- as.data.frame(t(test1.results))
  names(test1.results) <- names(comm_params)
  newnames <- paste(ids, seq(1:nrow(test1.results)), sep="")
  row.names(test1.results) <- newnames
  #test1.results <- round(test1.results)
  test1.results <- floor(test1.results) #gives better results
  results <- test1.results
}

simulate_comm.Norm <- function(comm_params, samples=100, ids="Norm_sim_comm", ...) {
  test1 <- lapply(seq_along(comm_params), function(i) {
    rnorm(samples, mean=comm_params[[i]][[1]][[1]],
          sd=comm_params[[i]][[1]][[2]])
  })
  test1.results <- do.call('rbind', test1)
  test1.results <- as.data.frame(t(test1.results))
  names(test1.results) <- names(comm_params)
  newnames <- paste(ids, seq(1:nrow(test1.results)), sep="")
  row.names(test1.results) <- newnames
  results <- test1.results
}

simulate_comm.CLR <- function(comm_params, samples=100, ids="CLR_sim_comm", ...) {
  test1 <- lapply(seq_along(comm_params), function(i) {
    rnorm(samples, mean=comm_params[[i]][[1]][[1]],
          sd=comm_params[[i]][[1]][[2]])
  })
  test1.results <- do.call('rbind', test1)
  test1.results <- as.data.frame(t(test1.results))
  names(test1.results) <- names(comm_params)
  newnames <- paste(ids, seq(1:nrow(test1.results)), sep="")
  row.names(test1.results) <- newnames
  results <- round(2^test1.results)
  return(results)
}

#######Need to debug below to have geomeans added automatically############
# simulate_comm.CLR <- function(comm_params, samples=100, ids="CLR_sim_comm", ...) {  #test version
#   test1 <- lapply(seq_along(comm_params), function(i) {
#     rnorm(samples, mean=comm_params[[i]][[1]][[1]],
#           sd=comm_params[[i]][[1]][[2]])
#   })
#   test1.results <- do.call('rbind', test1)
#   test1.results <- as.data.frame(t(test1.results))
#   names(test1.results) <- names(comm_params)
#   newnames <- paste(ids, seq(1:nrow(test1.results)), sep="")
#   row.names(test1.results) <- newnames
#   results <- round((2^(comm_params[[1]]$geomean+test1.results))-0.5) #use floor instead of round
#   return(results)
# }


