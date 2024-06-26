


Pj_TIRT_4_rank <- function(j,theta,aj=NULL,dj=NULL,BID=NULL,item.par=NULL){

  if(is.null(aj)){ #item.par and BID must be not NULL
    aj=item.par$a[BID$Block==j]
    dj=item.par$d[BID$Block==j]
  }

  if(is.vector(theta)){
    theta <- matrix(theta,nrow = 1)
  }
  if(!is.null(BID)){
    dim.j <- BID$Dim[BID$Block==j]
    thetaj = theta[,dim.j,drop=FALSE]
  }else{
    thetaj <- theta
  }

  alphatheta=aj*t(thetaj)
  p12=pnorm(alphatheta[1,]-alphatheta[2,]-dj[1]+dj[2])
  p13=pnorm(alphatheta[1,]-alphatheta[3,]-dj[1]+dj[3])
  p14=pnorm(alphatheta[1,]-alphatheta[4,]-dj[1]+dj[4])
  p21=1-p12
  p23=pnorm(alphatheta[2,]-alphatheta[3,]-dj[2]+dj[3])
  p24=pnorm(alphatheta[2,]-alphatheta[4,]-dj[2]+dj[4])
  p31=1-p13
  p32=1-p23
  p34=pnorm(alphatheta[3,]-alphatheta[4,]-dj[3]+dj[4])
  p41=1-p14
  p42=1-p24
  p43=1-p34
  pabcd=p12*p13*p14*p23*p24*p34
  pabdc=p12*p13*p14*p23*p24*p43
  pacbd=p12*p13*p14*p32*p24*p34
  pacdb=p12*p13*p14*p32*p42*p34
  padbc=p12*p13*p14*p23*p42*p43
  padcb=p12*p13*p14*p32*p42*p43

  pbacd=p21*p13*p14*p23*p24*p34
  pbadc=p21*p13*p14*p23*p24*p43
  pbcad=p21*p31*p14*p23*p24*p34
  pbcda=p21*p31*p41*p23*p24*p34
  pbdac=p21*p13*p41*p23*p24*p43
  pbdca=p21*p31*p41*p23*p24*p43

  pcabd=p12*p31*p14*p32*p24*p34
  pcadb=p12*p31*p14*p32*p42*p34
  pcbad=p21*p31*p14*p32*p24*p34
  pcbda=p21*p31*p41*p32*p24*p34
  pcdab=p12*p31*p41*p32*p42*p34
  pcdba=p21*p31*p41*p32*p42*p34

  pdabc=p12*p13*p41*p23*p42*p43
  pdacb=p12*p13*p41*p32*p42*p43
  pdbac=p21*p13*p41*p23*p42*p43
  pdbca=p21*p31*p41*p23*p42*p43
  pdcab=p12*p31*p41*p32*p42*p43
  pdcba=p21*p31*p41*p32*p42*p43
  prenorm=rbind(pabcd,pabdc,pacbd,pacdb,
                padbc,padcb,pbacd,pbadc,
                pbcad,pbcda,pbdac,pbdca,
                pcabd,pcadb,pcbad,pcbda,
                pcdab,pcdba,pdabc,pdacb,
                pdbac,pdbca,pdcab,pdcba)
  sumprenorm=matrix(colSums(prenorm)[col(prenorm)],nrow = 24)
  prenorm/sumprenorm
}

data.sim_TIRT_4_rank <- function(item.par,theta,BID){

  Y <- matrix(NA,nrow = nrow(theta),ncol = max(BID$Block))
  for(j in unique(BID$Block)){
    pj <- Pj_TIRT_4_rank(j=j,theta=theta,BID=BID,item.par=item.par)

    Y[,j] <- t(apply(pj,2,function(p)sample(x=nrow(pj),size=1,prob = p)))

  }
  Y
}




kernel_TIRT_4_rank <- function(Y,BID,a,d,sigm,theta,B=20,J,D,N,cor.matrix=FALSE,positive,fix.sigma=FALSE,hessian=FALSE){
# posterior for sampling theta - see Equation 4
arm.sampler_TIRT_4_rank <- function(thetai,p){
  logLi_TIRT_4_rank(thetai, a=a, d=d,yi=Y[i,], BID=BID,prior=TRUE,sigma=sigm)
}
# function to be optimized for estimating item parameters - see Equation 8
objfun_TIRT_4_rank <- function(x,thetaj,yj){
  -1*P.Yj_TIRT_4_rank(aj=(x[1:4]),dj=c(x[5:7],-sum(x[5:7])),thetaj,yj)
}
P.Yj_TIRT_4_rank <- function(aj,#vector of length 4
                            dj,#vector of length 4
                            thetaj,#matrix of N x 4
                            yj # vector of length N
){

  #alphatheta=matrix(rep(aj,each=nrow(thetaj)),ncol=nrow(thetaj),byrow=TRUE)*t(thetaj)
  alphatheta=aj*t(thetaj)
  p12=pnorm(alphatheta[1,]-alphatheta[2,]-dj[1]+dj[2])
  p13=pnorm(alphatheta[1,]-alphatheta[3,]-dj[1]+dj[3])
  p14=pnorm(alphatheta[1,]-alphatheta[4,]-dj[1]+dj[4])
  p21=1-p12
  p23=pnorm(alphatheta[2,]-alphatheta[3,]-dj[2]+dj[3])
  p24=pnorm(alphatheta[2,]-alphatheta[4,]-dj[2]+dj[4])
  p31=1-p13
  p32=1-p23
  p34=pnorm(alphatheta[3,]-alphatheta[4,]-dj[3]+dj[4])
  p41=1-p14
  p42=1-p24
  p43=1-p34
  pabcd=p12*p13*p14*p23*p24*p34
  pabdc=p12*p13*p14*p23*p24*p43
  pacbd=p12*p13*p14*p32*p24*p34
  pacdb=p12*p13*p14*p32*p42*p34
  padbc=p12*p13*p14*p23*p42*p43
  padcb=p12*p13*p14*p32*p42*p43

  pbacd=p21*p13*p14*p23*p24*p34
  pbadc=p21*p13*p14*p23*p24*p43
  pbcad=p21*p31*p14*p23*p24*p34
  pbcda=p21*p31*p41*p23*p24*p34
  pbdac=p21*p13*p41*p23*p24*p43
  pbdca=p21*p31*p41*p23*p24*p43

  pcabd=p12*p31*p14*p32*p24*p34
  pcadb=p12*p31*p14*p32*p42*p34
  pcbad=p21*p31*p14*p32*p24*p34
  pcbda=p21*p31*p41*p32*p24*p34
  pcdab=p12*p31*p41*p32*p42*p34
  pcdba=p21*p31*p41*p32*p42*p34

  pdabc=p12*p13*p41*p23*p42*p43
  pdacb=p12*p13*p41*p32*p42*p43
  pdbac=p21*p13*p41*p23*p42*p43
  pdbca=p21*p31*p41*p23*p42*p43
  pdcab=p12*p31*p41*p32*p42*p43
  pdcba=p21*p31*p41*p32*p42*p43
  prenorm=rbind(pabcd,pabdc,pacbd,pacdb,
                padbc,padcb,pbacd,pbadc,
                pbcad,pbcda,pbdac,pbdca,
                pcabd,pcadb,pcbad,pcbda,
                pcdab,pcdba,pdabc,pdacb,
                pdbac,pdbca,pdcab,pdcba)
  sumprenorm=matrix(colSums(prenorm)[col(prenorm)],nrow = 24)
  p.y=prenorm/sumprenorm
  p <- p.y[matrix(c(yj,seq_len(length(yj))),ncol = 2)]
  p[p<1e-5] <- 1e-5
  p[p>1-1e-5] <- 1-1e-5
  sum(log(p))
}

P.Yi_TIRT_4_rank <- function(a,#matrix 4 x J
                            d,#matrix 4 x J
                            theta, #matrix 4 x J  ability for each dimension on each item
                            yi# vector of J
){
  # yij is 1,...12
  alphatheta=a*theta
  p12=pnorm(alphatheta[1,]-alphatheta[2,]-d[1,]+d[2,])
  p13=pnorm(alphatheta[1,]-alphatheta[3,]-d[1,]+d[3,])
  p14=pnorm(alphatheta[1,]-alphatheta[4,]-d[1,]+d[4,])
  p21=1-p12
  p23=pnorm(alphatheta[2,]-alphatheta[3,]-d[2,]+d[3,])
  p24=pnorm(alphatheta[2,]-alphatheta[4,]-d[2,]+d[4,])
  p31=1-p13
  p32=1-p23
  p34=pnorm(alphatheta[3,]-alphatheta[4,]-d[3,]+d[4,])
  p41=1-p14
  p42=1-p24
  p43=1-p34
  pabcd=p12*p13*p14*p23*p24*p34
  pabdc=p12*p13*p14*p23*p24*p43
  pacbd=p12*p13*p14*p32*p24*p34
  pacdb=p12*p13*p14*p32*p42*p34
  padbc=p12*p13*p14*p23*p42*p43
  padcb=p12*p13*p14*p32*p42*p43

  pbacd=p21*p13*p14*p23*p24*p34
  pbadc=p21*p13*p14*p23*p24*p43
  pbcad=p21*p31*p14*p23*p24*p34
  pbcda=p21*p31*p41*p23*p24*p34
  pbdac=p21*p13*p41*p23*p24*p43
  pbdca=p21*p31*p41*p23*p24*p43

  pcabd=p12*p31*p14*p32*p24*p34
  pcadb=p12*p31*p14*p32*p42*p34
  pcbad=p21*p31*p14*p32*p24*p34
  pcbda=p21*p31*p41*p32*p24*p34
  pcdab=p12*p31*p41*p32*p42*p34
  pcdba=p21*p31*p41*p32*p42*p34

  pdabc=p12*p13*p41*p23*p42*p43
  pdacb=p12*p13*p41*p32*p42*p43
  pdbac=p21*p13*p41*p23*p42*p43
  pdbca=p21*p31*p41*p23*p42*p43
  pdcab=p12*p31*p41*p32*p42*p43
  pdcba=p21*p31*p41*p32*p42*p43
  prenorm=rbind(pabcd,pabdc,pacbd,pacdb,
                padbc,padcb,pbacd,pbadc,
                pbcad,pbcda,pbdac,pbdca,
                pcabd,pcadb,pcbad,pcbda,
                pcdab,pcdba,pdabc,pdacb,
                pdbac,pdbca,pdcab,pdcba)
  sumprenorm=matrix(colSums(prenorm)[col(prenorm)],nrow = 24)
  p.y=prenorm/sumprenorm
  p.y[matrix(c(yi,seq_len(length(yi))),ncol = 2)]
}

logLi_TIRT_4_rank <- function(thetai, #vector of length D  ability for person i
                             a,#matrix 4 x J
                             d,#matrix 4 x J
                             yi,# vector of J
                             BID,
                             prior=TRUE,
                             sigma=NULL){
  ll <- sum(log(P.Yi_TIRT_4_rank(a=a,d=d,theta=matrix(thetai[BID$Dim],nrow = 4),yi=yi)))
  lprior <- 0
  if(prior)
    lprior <- mvnfast::dmvn(X=matrix(thetai,nrow = 1),mu=rep(0,length(thetai)),sigma = sigma,log = TRUE)
  ll+lprior
}
  #number of parameters
  npar <- ifelse(fix.sigma,7*J,7*J+D*(D-1)/2)
  #matrix for all parameters x iterations
  pmatrix <- matrix(NA,npar,B)
  # theta.list <- list()
  for(b in seq_len(B)){
    # sampling theta - Equation 4
    theta <- foreach (i=seq_len(N),.combine="rbind",.errorhandling = "stop",.inorder = TRUE)%dopar%{

                        armspp::arms_gibbs(n_samples = 1,previous = theta[i,],log_pdf = arm.sampler_TIRT_4_rank,
                                           lower = rep(-6,D),upper = rep(6,D),metropolis = TRUE)
                      }
    if(!fix.sigma){
      # if estimating sigma
      if(cor.matrix){
        # use correlation matrix as the estimate if cor.matrix = TRUE
        sigm <- cor(theta)
      }else{
        # estimate sigma using the proximal gradient descent algorithm - see Equations 12 - 15
        sigm <- try(lvmcomp:::calcu_sigma_cmle_cpp(theta,1e-5),silent=TRUE)
        if(class(sigm)[1]=="try-error"){
          sigm <- cor(theta)
          cat("\nEstimated covariance matrix is not positive definite.")
        }
      }
    }


    # estimating item parameters using L-BFGS-B algorithm: see Equation 8
    opt <- foreach (j=1:J,.combine="cbind",.errorhandling = "stop",.inorder = TRUE)%dopar%{
                      lo <- c(ifelse(positive[which(BID$Block==j)]>0,c(0,0,0,0),c(-Inf,-Inf,-Inf,-Inf)),-6,-6,-6)
                      up <- c(ifelse(positive[which(BID$Block==j)]>0,c(Inf,Inf,Inf,Inf),c(0,0,0,0)),6,6,6)
                      optim(par=c(a[,j],d[1:3,j]),fn = objfun_TIRT_4_rank,method = "L-BFGS-B",lower = lo, upper = up,
                            thetaj=theta[,BID$Dim[which(BID$Block==j)]],yj=Y[,j],control = list(maxit=20))$par
                    }
    hess=diag(rep(0,J*7))
    if(hessian==TRUE){
      hess_vec=foreach (j=1:J,.combine="cbind",.errorhandling = "stop",.inorder = TRUE)%dopar%{
        lo <- c(ifelse(positive[which(BID$Block==j)]>0,c(0,0,0,0),c(-Inf,-Inf,-Inf,-Inf)),-6,-6,-6)
        up <- c(ifelse(positive[which(BID$Block==j)]>0,c(Inf,Inf,Inf,Inf),c(0,0,0,0)),6,6,6)
        optim(par=c(a[,j],d[1:3,j]),fn = objfun_TIRT_4_rank,method = "L-BFGS-B",lower = lo, upper = up,
              thetaj=theta[,BID$Dim[which(BID$Block==j)]],yj=Y[,j],control = list(maxit=20),hessian=TRUE)$hessian
      }
      for (i in seq(1, J*7, by = 7)) {
        start_row <- i
        start_col <- i
        end_row <- start_row + 6
        end_col <- start_col + 6

        hess[start_row:end_row, start_col:end_col] <- hess_vec[, ((start_row - 1) / 7 * 7 + 1):((start_row - 1) / 7 * 7 + 7)]
      }
    }



    a <- opt[1:4,]
    d <- rbind(opt[5:7,],-1*colSums(opt[5:7,]))
    if(fix.sigma){
      pmatrix[,b] <- c(opt)
    }else{
      pmatrix[,b] <- c(c(opt),c(sigm[lower.tri(sigm)]))
    }


  }

  return(list(pmatrix=pmatrix,a=a,d=d,sigm=sigm,theta=theta,hess=hess))
}



param.vec.2.ads_TIRT_4_rank <- function(pv,J,D,sigm=NULL,fix.sigma=FALSE){
  #pv: the vector of all parameters (a,d,sig)
  ad <- matrix(pv[1:(7*J)],nrow = 7)
  a <- ad[1:4,]
  d <- rbind(ad[5:7,],-1*colSums(ad[5:7,]))
  if(!fix.sigma){
    sig <- diag(D)
    sig[lower.tri(sig)] <- pv[((7*J)+1):length(pv)]
    sigm=Matrix::forceSymmetric(sig,uplo="L")
  }

  list(a=a,d=d,sigm=sigm)
}
logLi_TIRT_4_rank <- function(thetai, #vector of length D  ability for person i
                              a,#matrix 4 x J
                              d,#matrix 4 x J
                              yi,# vector of J
                              BID,
                              prior=TRUE,
                              sigma=NULL){
  ll <- sum(log(P.Yi_TIRT_4_rank(a=a,d=d,theta=matrix(thetai[BID$Dim],nrow = 4),yi=yi)))
  lprior <- 0
  if(prior)
    lprior <- mvnfast::dmvn(X=matrix(thetai,nrow = 1),mu=rep(0,length(thetai)),sigma = sigma,log = TRUE)
  ll+lprior
}

P.Yi_TIRT_4_rank <- function(a,#matrix 4 x J
                             d,#matrix 4 x J
                             theta, #matrix 4 x J  ability for each dimension on each item
                             yi# vector of J
){
  # yij is 1,...12
  alphatheta=a*theta
  p12=pnorm(alphatheta[1,]-alphatheta[2,]-d[1,]+d[2,])
  p13=pnorm(alphatheta[1,]-alphatheta[3,]-d[1,]+d[3,])
  p14=pnorm(alphatheta[1,]-alphatheta[4,]-d[1,]+d[4,])
  p21=1-p12
  p23=pnorm(alphatheta[2,]-alphatheta[3,]-d[2,]+d[3,])
  p24=pnorm(alphatheta[2,]-alphatheta[4,]-d[2,]+d[4,])
  p31=1-p13
  p32=1-p23
  p34=pnorm(alphatheta[3,]-alphatheta[4,]-d[3,]+d[4,])
  p41=1-p14
  p42=1-p24
  p43=1-p34
  pabcd=p12*p13*p14*p23*p24*p34
  pabdc=p12*p13*p14*p23*p24*p43
  pacbd=p12*p13*p14*p32*p24*p34
  pacdb=p12*p13*p14*p32*p42*p34
  padbc=p12*p13*p14*p23*p42*p43
  padcb=p12*p13*p14*p32*p42*p43

  pbacd=p21*p13*p14*p23*p24*p34
  pbadc=p21*p13*p14*p23*p24*p43
  pbcad=p21*p31*p14*p23*p24*p34
  pbcda=p21*p31*p41*p23*p24*p34
  pbdac=p21*p13*p41*p23*p24*p43
  pbdca=p21*p31*p41*p23*p24*p43

  pcabd=p12*p31*p14*p32*p24*p34
  pcadb=p12*p31*p14*p32*p42*p34
  pcbad=p21*p31*p14*p32*p24*p34
  pcbda=p21*p31*p41*p32*p24*p34
  pcdab=p12*p31*p41*p32*p42*p34
  pcdba=p21*p31*p41*p32*p42*p34

  pdabc=p12*p13*p41*p23*p42*p43
  pdacb=p12*p13*p41*p32*p42*p43
  pdbac=p21*p13*p41*p23*p42*p43
  pdbca=p21*p31*p41*p23*p42*p43
  pdcab=p12*p31*p41*p32*p42*p43
  pdcba=p21*p31*p41*p32*p42*p43
  prenorm=rbind(pabcd,pabdc,pacbd,pacdb,
                padbc,padcb,pbacd,pbadc,
                pbcad,pbcda,pbdac,pbdca,
                pcabd,pcadb,pcbad,pcbda,
                pcdab,pcdba,pdabc,pdacb,
                pdbac,pdbca,pdcab,pdcba)
  sumprenorm=matrix(colSums(prenorm)[col(prenorm)],nrow = 24)
  p.y=prenorm/sumprenorm
  p.y[matrix(c(yi,seq_len(length(yi))),ncol = 2)]
}
logL_TIRT_4_rank <- function(theta, #vector of length D  ability for person i
                             a,#matrix 3 x J
                             d,#matrix 3 x J
                             yi,# vector of J
                             BID,
                             sigma=NULL){
  # yij is 1,...6
  thetai <- matrix(theta[BID$Dim],nrow = 4)
  alphatheta=a*thetai
  p12=pnorm(alphatheta[1,]-alphatheta[2,]-d[1,]+d[2,])
  p13=pnorm(alphatheta[1,]-alphatheta[3,]-d[1,]+d[3,])
  p14=pnorm(alphatheta[1,]-alphatheta[4,]-d[1,]+d[4,])
  p21=1-p12
  p23=pnorm(alphatheta[2,]-alphatheta[3,]-d[2,]+d[3,])
  p24=pnorm(alphatheta[2,]-alphatheta[4,]-d[2,]+d[4,])
  p31=1-p13
  p32=1-p23
  p34=pnorm(alphatheta[3,]-alphatheta[4,]-d[3,]+d[4,])
  p41=1-p14
  p42=1-p24
  p43=1-p34
  pabcd=p12*p13*p14*p23*p24*p34
  pabdc=p12*p13*p14*p23*p24*p43
  pacbd=p12*p13*p14*p32*p24*p34
  pacdb=p12*p13*p14*p32*p42*p34
  padbc=p12*p13*p14*p23*p42*p43
  padcb=p12*p13*p14*p32*p42*p43

  pbacd=p21*p13*p14*p23*p24*p34
  pbadc=p21*p13*p14*p23*p24*p43
  pbcad=p21*p31*p14*p23*p24*p34
  pbcda=p21*p31*p41*p23*p24*p34
  pbdac=p21*p13*p41*p23*p24*p43
  pbdca=p21*p31*p41*p23*p24*p43

  pcabd=p12*p31*p14*p32*p24*p34
  pcadb=p12*p31*p14*p32*p42*p34
  pcbad=p21*p31*p14*p32*p24*p34
  pcbda=p21*p31*p41*p32*p24*p34
  pcdab=p12*p31*p41*p32*p42*p34
  pcdba=p21*p31*p41*p32*p42*p34

  pdabc=p12*p13*p41*p23*p42*p43
  pdacb=p12*p13*p41*p32*p42*p43
  pdbac=p21*p13*p41*p23*p42*p43
  pdbca=p21*p31*p41*p23*p42*p43
  pdcab=p12*p31*p41*p32*p42*p43
  pdcba=p21*p31*p41*p32*p42*p43
  prenorm=rbind(pabcd,pabdc,pacbd,pacdb,
                padbc,padcb,pbacd,pbadc,
                pbcad,pbcda,pbdac,pbdca,
                pcabd,pcadb,pcbad,pcbda,
                pcdab,pcdba,pdabc,pdacb,
                pdbac,pdbca,pdcab,pdcba)
  sumprenorm=matrix(colSums(prenorm)[col(prenorm)],nrow = 24)
  p.y=(prenorm/sumprenorm)[matrix(c(yi,seq_len(length(yi))),ncol = 2)]
  if(is.null(sigma)){
    # log likelihood
    logP <- sum(log(p.y))
  }else{
    # log likelihood + log prior  dmvn-多元概率密度函数：根据先验分布找到当前能力组合的概率值，取对数
    # 关于dmvn函数，该函数底层代码为C++，请参考“先验信息相关函数文件夹”内的dmvtCpp，与之相关的内容均在本文件夹下
    logP <- sum(log(p.y)) + dmvn(theta,mu=rep(0,length(theta)),sigma = sigma,log = TRUE)
  }
  -1*logP
}
