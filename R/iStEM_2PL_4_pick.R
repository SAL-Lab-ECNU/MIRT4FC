

iStEM_2PL_4_pick <- function(Y,BID,positive=rep(TRUE,nrow(BID)),M=10,B=20,SE='Bootstrap',blocksize=4,
                 sigma=NULL,theta=NULL,fix.sigma = FALSE,burnin.maxitr=40,model='2PL',
                 maxitr=500,eps1=1.5,eps2=0.4,frac1=.2,frac2=.5,cores=NULL,h=NULL){

  ##################################
  #
  #
  #
  #
  ##################################


  s1 <- Sys.time()
  #number of dimensions
  D <- length(unique(BID$Dim))
  N <- nrow(Y) #number of participants
  J <- max(BID$Block) # number of blocks
  Q <- nrow(BID) # number of questions
  IN <- 2*blocksize-1 # Number of item parameter
  npar <- IN*J+D*(D-1)/2
  #initial item parameters
  a <- matrix(positive,nrow=blocksize,ncol=J)
  d <- matrix(rnorm(blocksize*J),nrow=blocksize)
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
  pplist <- list()
  hess=list()
  ###########################
  #
  # burn-in phase
  #
  ###########################
  total.number.of.batch <- burn.in.size <- 0
  # the initial MxB iterations
  cat("\nBurn-in phase:")
  # the kernel_4_pick function returns item parameter estimates, person parameter estimates and sigma estimates
  x <- kernel_2PL_4_pick(Y=Y,BID=BID,a=a,d=d,sigm=sigm,theta=theta,B=B,J=J,D=D,N=N,cor.matrix = TRUE,positive=positive,fix.sigma=fix.sigma,hessian = TRUE)
  a <- x$a
  d <- x$d
  sigm <- x$sigm
  theta <- x$theta
  for(m in 1:M){

    cat("\n  # of batch = ",m)
    x <- kernel_2PL_4_pick(Y=Y,BID=BID,a=a,d=d,sigm=sigm,theta=theta,B=B,J=J,D=D,N=N,cor.matrix = (m<M/2),positive=positive,fix.sigma=fix.sigma,hessian = TRUE)
    plist[[m]] <- x$pmatrix
    pplist <- c(pplist,x$ptheta)
    hess=c(hess,x$hess_list)
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
  #print(z)
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

    x <- kernel_2PL_4_pick(Y=Y,BID=BID,a=a,d=d,sigm=sigm,theta=theta,B=B,J=J,D=D,N=N,positive=positive,fix.sigma=fix.sigma,hessian = TRUE)
    a <- x$a
    d <- x$d
    sigm <- x$sigm
    #print(sigm)
    theta <- x$theta
    plist <- plist[-1]
    plist[[M]] <- x$pmatrix
    pplist <- pplist[-(1:B)]
    pplist <- c(pplist,x$ptheta)
    hess=hess[-(1:B)]
    hess=c(hess,x$hess_list)
    d.hat <- batch.var(plist,n=M)
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
  # m <- M
  while(sum(z^2)>=npar*eps1&&m<burnin.maxitr&&max(d.hat*N)>=eps2){ # if the chain is not stable....
    if(m==M){
      cat("  sum z^2 / npar = ",sum(z^2)/npar, " criterion = ",eps1," max delta hat = ",max(d.hat)*N, " criterion = ",eps2)
    }else{
      cat("\n  # of batch = ",m," sum z^2 / npar = ",sum(z^2)/npar, " criterion = ",eps1," max delta hat = ",max(d.hat)*N, " criterion = ",eps2)
    }

    x <- kernel_2PL_4_pick(Y=Y,BID=BID,a=a,d=d,sigm=sigm,theta=theta,B=B,J=J,D=D,N=N,positive=positive,fix.sigma=fix.sigma,hessian = TRUE)
    a <- x$a
    d <- x$d
    sigm <- x$sigm
    #print(sigm)
    theta <- x$theta
    plist <- plist[-1]
    plist[[M]] <- x$pmatrix
    pplist <- pplist[-(1:B)]
    pplist <- c(pplist,x$ptheta)
    hess=hess[-(1:B)]
    hess=c(hess,x$hess_list)
    d.hat <- batch.var(plist,n=M)
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

    x <- kernel_2PL_4_pick(Y=Y,BID=BID,a=a,d=d,sigm=sigm,theta=theta,B=B,J=J,D=D,N=N,positive=positive,fix.sigma=fix.sigma,hessian = TRUE)
    a <- x$a
    d <- x$d
    sigm <- x$sigm
    theta <- x$theta
    plist[[n]] <- x$pmatrix
    pplist <- c(pplist,x$ptheta)
    hess=c(hess,x$hess_list)
    d.hat <- batch.var(plist,n=n)
  }
    cat("\n  # of valid batch = ",n," max delta hat = ",max(d.hat)*N, " criterion = ",eps2)

  cat("\nLenth of final MC chain = ",n*B,"\n")
  total.number.of.batch <- total.number.of.batch + n - M
  est <- param.vec.2.ads_2PL_4_pick(rowMeans(Reduce(cbind,plist)),J,D,sigm=sigm,fix.sigma=fix.sigma)
  matrix_array <- array(unlist(hess), dim = c(ncol(hess[[1]]), ncol(hess[[1]]), length(hess)))
  hess_mean <- apply(matrix_array, c(1, 2), mean)
  se.item=list(a=est$a,d=est$d,itemlist=plist,thetalist=pplist)
  #计算SE

  if(SE=='Bootstrap'){
    SE_item=itemSE_bootstrap(plist=plist,blocksize=blocksize,IN=IN,model=model,J=J,D=D)
  }else if(SE=='FDM'|SE=='XPD'|SE=='CDM'|SE=='RES'|SE=='Louis'|SE=='Sandwich'){
    SE_item=itemSE_Simple_2PL_4_pick(se.item,BID,Y,SE,IN,hess_mean,h1=h)
    SE_item$SE.a=t(SE_item$SE.a)
    SE_item$SE.d=t(SE_item$SE.d)
  }else if(SE=='complete'){
    se <- matrix(sqrt(diag(solve(hess_mean))),nrow=max(BID$Block),ncol=7,byrow = TRUE)
    se.a <- se[,1:4]
    se.d <- se[,5:7]
    se.d <- cbind(se.d,sqrt(se.d[,1]^2+se.d[,2]^2+se.d[,3]^2))
    SE_item=list(SE.a=se.a,SE.d=se.d,SE.sigm=NULL)
    SE_item$SE.a=t(SE_item$SE.a)
    SE_item$SE.d=t(SE_item$SE.d)
  }
  stopCluster(cl)
  s2 <- Sys.time()
  return(list(SE.type=SE,a=est$a,SE.a=SE_item$SE.a,d=est$d,SE.d=SE_item$SE.d,
              sigm=est$sigm,SE.sigm=SE_item$SE.sigm,total.number.of.batch=total.number.of.batch,
              final.chain=n*B,burn.in.size=burn.in.size,itemlist=plist,thetalist=pplist,
              timeused = s2 - s1, start.time = s1, end.time = s2,hess=hess))
}




