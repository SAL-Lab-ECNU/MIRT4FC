#' stEM_4_pick
#'
#' @param Y          A # of subjects x # of blocks matrix; item responses
#' @param BID        A # of statements x 3 matrix; item information,
#'     columns are "Block", "Item" and "Dimensions"
#' @param positive   A logical vector; indicating whether each statement is
#'     positive directional or not
#' @param M          A number; # of batch
#' @param B          A number; # of iterations in each batch
#' @param a          A vector; length = # of statements, initial alpha parameters
#' @param d          A vector; length = # of statements, initial beta parameters
#' @param item.par   A data frame; initial parameters for a and d
#' @param sigma      A # of dimensions x # of dimensions matrix; initial sigma
#'     parameters
#' @param theta      A # of subjects x # of dimensions matrix; initial theta parameters
#' @param fix.sigma  Logical; TRUE if sigma is estimated
#' @param burnin.maxitr   A number; max burn-in allowed
#' @param maxitr     A number; max iterations allowed
#' @param eps1       A number; stability criteria
#' @param eps2       A number; convergence criterion
#' @param frac1      A number; cutoffs for calculating Geweke z
#' @param frac2      A number; cutoffs for calculating Geweke z
#' @param cores      A number; number of cores
#'
#' @return A list of data; including estimated a, d, and sigm parameters,
#'    total batch number, final chain size, burn-in size, time
#'
#'

stEM_4_pick <- function(Y,BID,positive=rep(TRUE,nrow(BID)),M=10,B=20,a=NULL,d=NULL,
                 item.par=NULL,sigma=NULL,theta=NULL,fix.sigma = FALSE,burnin.maxitr=40,
                 maxitr=500,eps1=1.5,eps2=0.4,frac1=.2,frac2=.5,cores=NULL){

  ##################################
  #
  #
  #
  #
  ##################################

  required.packages <- c("armspp", "doParallel","foreach","doRNG","coda","mvnfast","lvmcomp")
  to.be.installed <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  if(length(to.be.installed)) install.packages(to.be.installed)
  stopifnot(sapply(required.packages,require,character.only=TRUE))
  s1 <- Sys.time()
  #number of dimensions
  D <- length(unique(BID$Dim))
  N <- nrow(Y) #number of participants
  J <- max(BID$Block) # number of blocks
  Q <- nrow(BID) # number of questions
  npar <- 7*J+D*(D-1)/2
  #initial a parameters
  if(is.null(a)){
    if(!is.null(item.par)){
      a <- matrix(item.par$a,nrow = 4) #4 x J
    }else{
      a <- matrix(positive,nrow=4,ncol=J)
    }
  }
  #initial d parameters
  if(is.null(d)){
    if(!is.null(item.par)){
      d <- matrix(item.par$d,nrow = 4) #4 x J
    }else{
      d <- matrix(rnorm(4*J),nrow=4)
    }

  }
  #initial theta
  if(is.null(theta)){
    theta <- matrix(0,nrow = N,ncol = D)
  }
  #initial covariance matrix
  if(is.null(sigma)){
    sigm <- diag(D)
  }else{
    sigm <- sigma
  }

  # parallel computing settings

  nCPUcores = detectCores()
  if (nCPUcores < 3) {
    registerDoSEQ()
  }else{
    if(is.null(cores)){
      cl = makeCluster(nCPUcores-1)
      registerDoParallel(cl)
    }else if(cores<nCPUcores){
      cl = makeCluster(cores)
      registerDoParallel(cl)
    }else{
      cl = makeCluster(nCPUcores-1)
      registerDoParallel(cl)
    }
  }
  plist <- list()

  ###########################
  #
  # burn-in phase
  #
  ###########################
  total.number.of.batch <- burn.in.size <- 0
  # the initial MxB iterations
  cat("\nBurn-in phase:")
  # the kernel_4_pick function returns item parameter estimates, person parameter estimates and sigma estimates
  x <- kernel_4_pick(Y=Y,BID=BID,a=a,d=d,sigm=sigm,theta=theta,B=B,J=J,D=D,N=N,cor.matrix = TRUE,positive=positive,fix.sigma=fix.sigma)
  a <- x$a
  d <- x$d
  sigm <- x$sigm
  theta <- x$theta
  for(m in 1:M){

    cat("\n  # of batch = ",m)
    x <- kernel_4_pick(Y=Y,BID=BID,a=a,d=d,sigm=sigm,theta=theta,B=B,J=J,D=D,N=N,cor.matrix = (m<M/2),positive=positive,fix.sigma=fix.sigma)
    plist[[m]] <- x$pmatrix
    a <- x$a
    d <- x$d
    sigm <- x$sigm
    theta <- x$theta
    # print(sigm)
    if(m>=burnin.maxitr) break
  }
  #terminates burn-in using Geweke statistic z? --- see Page 3 of the manual
  mcmc.par <- coda::mcmc(t(Reduce(cbind,plist)))
  z <- coda::geweke.diag(mcmc.par,frac1 = frac1,frac2 = frac2)$z
  print(z)
  for(i in 1:npar){
    if(is.nan(z[i])){
      z[i]=0
    }
  }
d.hat <- batch.var(plist,n=m)
 for(i in 1:5){ # 5 iterations
    if(m==M){
      cat("  sum z^2 / npar = ",sum(z^2)/npar, " criterion = ",eps1," max delta hat = ",max(d.hat)*N, " criterion = ",eps2)
    }else{
      cat("\n  # of batch = ",m," sum z^2 / npar = ",sum(z^2)/npar, " criterion = ",eps1," max delta hat = ",max(d.hat)*N, " criterion = ",eps2)
    }

    x <- kernel_4_pick(Y=Y,BID=BID,a=a,d=d,sigm=sigm,theta=theta,B=B,J=J,D=D,N=N,positive=positive,fix.sigma=fix.sigma)
    a <- x$a
    d <- x$d
    sigm <- x$sigm
    #print(sigm)
    theta <- x$theta
    plist <- plist[-1]
    plist[[M]] <- x$pmatrix
    d.hat <- batch.var(plist,n=m)
    #terminates burn-in using Geweke statistic z? --- see Page 3 of the manual
    mcmc.par <- coda::mcmc(t(Reduce(cbind,plist)))
    z <- coda::geweke.diag(mcmc.par,frac1 = frac1,frac2 = frac2)$z
    cat("  sum z^2 / npar = ",sum(z^2)/npar)
    burn.in.size <- burn.in.size + B
    m <- m+1
  }
  # m <- M
  while(sum(z^2)>=npar*eps1&&m<burnin.maxitr&&max(d.hat*N)>=eps2){ # if the chain is not stable....
    if(m==M){
      cat("  sum z^2 / npar = ",sum(z^2)/npar, " criterion = ",eps1," max delta hat = ",max(d.hat)*N, " criterion = ",eps2)
    }else{
      cat("\n  # of batch = ",m," sum z^2 / npar = ",sum(z^2)/npar, " criterion = ",eps1," max delta hat = ",max(d.hat)*N, " criterion = ",eps2)
    }

    x <- kernel_4_pick(Y=Y,BID=BID,a=a,d=d,sigm=sigm,theta=theta,B=B,J=J,D=D,N=N,positive=positive,fix.sigma=fix.sigma)
    a <- x$a
    d <- x$d
    sigm <- x$sigm
    #print(sigm)
    theta <- x$theta
    plist <- plist[-1]
    plist[[M]] <- x$pmatrix
    d.hat <- batch.var(plist,n=m)
    #terminates burn-in using Geweke statistic z? --- see Page 3 of the manual
    mcmc.par <- coda::mcmc(t(Reduce(cbind,plist)))
    z <- coda::geweke.diag(mcmc.par,frac1 = frac1,frac2 = frac2)$z
    for(i in 1:npar){
      if(is.nan(z[i])){
        z[i]=0
      }
    }
    cat("  sum z^2 / npar = ",sum(z^2)/npar)
    burn.in.size <- burn.in.size + B
    m <- m+1
  }

  cat("\n  # of batch = ",m," sum z^2 / npar = ",sum(z^2)/npar, " criterion = ",eps1," max delta hat = ",max(d.hat)*N, " criterion = ",eps2)
  total.number.of.batch <- m
  cat("\nBurn in batch = ",1 + burn.in.size/B," burn-in iterations = ",B + burn.in.size)

  ###########################
  #
  # determining chain length
  #
  ###########################
  n <- M
  d.hat <- batch.var(plist,n=n)
    cat("\nAfter burn-in phase:")
  while(max(d.hat*N)>=eps2&&n<maxitr){
    cat("\n  # of valid batch = ",n," max delta hat = ",max(d.hat)*N, " criterion = ",eps2)
    n <- n+1

    x <- kernel_4_pick(Y=Y,BID=BID,a=a,d=d,sigm=sigm,theta=theta,B=B,J=J,D=D,N=N,positive=positive,fix.sigma=fix.sigma)
    a <- x$a
    d <- x$d
    sigm <- x$sigm
    theta <- x$theta
    plist[[n]] <- x$pmatrix

    d.hat <- batch.var(plist,n=n)
  }
    cat("\n  # of valid batch = ",n," max delta hat = ",max(d.hat)*N, " criterion = ",eps2)

  cat("\nLenth of final MC chain = ",n*B)
  total.number.of.batch <- total.number.of.batch + n - M
  est <- param.vec.2.ads_4_pick(rowMeans(Reduce(cbind,plist)),J,D,sigm=sigm,fix.sigma=fix.sigma)

  stopCluster(cl)
  s2 <- Sys.time()
  return(list(a=est$a,d=est$d,sigm=est$sigm,total.number.of.batch=total.number.of.batch,
              final.chain=n*B,burn.in.size=burn.in.size,plist=plist,
              timeused = s2 - s1, start.time = s1, end.time = s2))
}


