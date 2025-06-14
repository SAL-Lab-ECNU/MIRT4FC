


Pj_2PL_3_rank <- function(j,theta,aj=NULL,dj=NULL,BID=NULL,item.par=NULL){
  response.pattern <- matrix(c(1,2,3,
                               1,3,2,
                               2,1,3,
                               2,3,1,
                               3,1,2,
                               3,2,1),ncol = 3,byrow = TRUE)
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


  z <- exp(aj*t(thetaj)-dj) # 3 x N
  sumz <- matrix(colSums(z)[col(z)],nrow = 3)
  red.sumz <- sumz-z
  pj.a.abc <- z/sumz #three rows: za/sum(z),zb/sum(z),zc/sum(z); J columns
  pj.b.bc <- z[response.pattern[,2],,drop=FALSE]/red.sumz[response.pattern[,1],,drop=FALSE]
  pj.a.abc[response.pattern[,1],,drop=FALSE]*pj.b.bc
}

data.sim_2PL_3_rank <- function(item.par,theta,BID){

  Y <- matrix(NA,nrow = nrow(theta),ncol = max(BID$Block))
  for(j in unique(BID$Block)){
    pj <- Pj_2PL_3_rank(j=j,theta=theta,BID=BID,item.par=item.par)

    Y[,j] <- t(apply(pj,2,function(p)sample(x=nrow(pj),size=1,prob = p)))

  }
  Y
}


kernel_2PL_3_rank <- function(Y,BID,a,d,sigm,theta,B=20,J,D,N,cor.matrix=FALSE,positive,fix.sigma=FALSE,hessian=TRUE){
# posterior for sampling theta - see Equation 4
arm.sampler_2PL_3_rank <- function(thetai,p){
  logLi_2PL_3_rank(thetai, a=a, d=d,yi=Y[i,], BID=BID,prior=TRUE,sigma=sigm)
}
# function to be optimized for estimating item parameters - see Equation 8
objfun_2PL_3_rank <- function(x,thetaj,yj){
  -1*P.Yj_2PL_3_rank(aj=(x[1:3]),dj=c(x[4:5],-sum(x[4:5])),thetaj,yj)
}

logLi_2PL_3_rank <- function(thetai, #vector of length D  ability for person i
                             a,#matrix 3 x J
                             d,#matrix 3 x J
                             yi,# vector of J
                             BID,
                             prior=TRUE,
                             sigma=NULL){
  ll <- sum(log(P.Yi_2PL_3_rank(a=a,d=d,theta=matrix(thetai[BID$Dim],nrow = 3),yi=yi)))
  lprior <- 0
  if(prior)
    lprior <- mvnfast::dmvn(X=matrix(thetai,nrow = 1),mu=rep(0,length(thetai)),sigma = sigma,log = TRUE)
  ll+lprior
}


P.Yj_2PL_3_rank <- function(aj,#vector of length 3
                            dj,#vector of length 3
                            thetaj,#matrix of N x 3
                            yj # vector of length N
){

  response.pattern <- matrix(c(1,2,3,
                               1,3,2,
                               2,1,3,
                               2,3,1,
                               3,1,2,
                               3,2,1),ncol = 3,byrow = TRUE)
  z <- exp(aj*t(thetaj)-dj) # 3 x N
  sumz <- matrix(colSums(z)[col(z)],nrow = 3)
  red.sumz <- sumz-z
  p.a.abc <- z/sumz #three rows: za/sum(z),zb/sum(z),zc/sum(z); J columns
  p.b.bc <- z[response.pattern[,2],]/red.sumz[response.pattern[,1],]
  p.y <- p.a.abc[response.pattern[,1],]*p.b.bc
  p <- p.y[matrix(c(yj,seq_len(length(yj))),ncol = 2)]
  p[p<1e-5] <- 1e-5
  p[p>1-1e-5] <- 1-1e-5
  sum(log(p))
}

P.Yi_2PL_3_rank <- function(a,#matrix 3 x J
                            d,#matrix 3 x J
                            theta, #matrix 3 x J  ability for each dimension on each item
                            yi# vector of J
){
  # yij is 1,...6
  response.pattern <- matrix(c(1,2,3,
                               1,3,2,
                               2,1,3,
                               2,3,1,
                               3,1,2,
                               3,2,1),ncol = 3,byrow = TRUE)
  z <- exp(a*theta-d) #3 x J
  sumz <- matrix(colSums(z)[col(z)],nrow = 3)
  red.sumz <- sumz-z
  p.a.abc <- z/sumz #three rows: za/sum(z),zb/sum(z),zc/sum(z); J columns
  p.b.bc <- z[response.pattern[,2],]/red.sumz[response.pattern[,1],]
  p.y <- p.a.abc[response.pattern[,1],]*p.b.bc

  p.y[matrix(c(yi,seq_len(length(yi))),ncol = 2)]
}
  #number of parameters
  npar <- ifelse(fix.sigma,5*J,5*J+D*(D-1)/2)
  #matrix for all parameters x iterations
  pmatrix <- matrix(NA,npar,B)
  # theta.list <- list()
  ptheta=list()
  hess_list=list()
  for(b in seq_len(B)){
    # sampling theta - Equation 4
    theta <- foreach (i=seq_len(N),.combine="rbind",.errorhandling = "stop",.inorder = TRUE)%dopar%{

                        armspp::arms_gibbs(n_samples = 1,previous = theta[i,],log_pdf = arm.sampler_2PL_3_rank,
                                           lower = rep(-6,D),upper = rep(6,D),metropolis = fix.sigma)
                      }
    if (!fix.sigma) {
      if (cor.matrix) {
        sigm <- cor(theta)
      } else {
        # Attempt to estimate sigma using lvmcomp:::calcu_sigma_cmle_cpp()
        sigm <- tryCatch(
          {
            lvmcomp:::calcu_sigma_cmle_cpp(theta, 1e-5)
          },
          error = function(e) {
            cat("\nEstimated covariance matrix is not positive definite,fall back to using correlation matrix")
            #
            cor(theta)
          }
        )
      }
    }
    # estimating item parameters using L-BFGS-B algorithm: see Equation 8

    opt <- foreach (j=1:J,.combine="cbind",.errorhandling = "stop",.inorder = TRUE)%dopar%{
                      lo <- c(ifelse(positive[which(BID$Block==j)]>0,c(0,0,0),c(-Inf,-Inf,-Inf)),-6,-6)
                      up <- c(ifelse(positive[which(BID$Block==j)]>0,c(Inf,Inf,Inf),c(0,0,0)),6,6)
                      optim(par=c(a[,j],d[1:2,j]),fn = objfun_2PL_3_rank,method = "L-BFGS-B",lower = lo, upper = up,
                            thetaj=theta[,BID$Dim[which(BID$Block==j)]],yj=Y[,j],control = list(maxit=20))$par
    }
    hess=diag(rep(0,J*5))
    if(hessian==TRUE){
      hess_vec=foreach (j=1:J,.combine="cbind",.errorhandling = "stop",.inorder = TRUE)%dopar%{
        lo <- c(ifelse(positive[which(BID$Block==j)]>0,c(0,0,0),c(-Inf,-Inf,-Inf)),-6,-6)
        up <- c(ifelse(positive[which(BID$Block==j)]>0,c(Inf,Inf,Inf),c(0,0,0)),6,6)
        optim(par=c(a[,j],d[1:2,j]),fn = objfun_2PL_3_rank,method = "L-BFGS-B",lower = lo, upper = up,
              thetaj=theta[,BID$Dim[which(BID$Block==j)]],yj=Y[,j],control = list(maxit=20),hessian=TRUE)$hessian
      }
      for (i in seq(1, J*5, by = 5)) {
        start_row <- i
        start_col <- i
        end_row <- start_row + 4
        end_col <- start_col + 4

        hess[start_row:end_row, start_col:end_col] <- hess_vec[, ((start_row - 1) / 5 * 5 + 1):((start_row - 1) / 5 * 5 + 5)]
      }
    }


    a <- opt[1:3,]
    d <- rbind(opt[4:5,],-1*colSums(opt[4:5,]))
    if(fix.sigma){
      pmatrix[,b] <- c(opt)
    }else{
      pmatrix[,b] <- c(c(opt),c(sigm[lower.tri(sigm)]))
    }
    ptheta[[b]]=theta
    hess_list[[b]]=hess
  }

  return(list(pmatrix=pmatrix,a=a,d=d,sigm=sigm,theta=theta,hess_list=hess_list,ptheta=ptheta))
}



param.vec.2.ads_2PL_3_rank <- function(pv,J,D,sigm=NULL,fix.sigma=FALSE){
  #pv: the vector of all parameters (a,d,sig)
  ad <- matrix(pv[1:(5*J)],nrow = 5)
  a <- ad[1:3,]
  d <- rbind(ad[4:5,],-1*colSums(ad[4:5,]))
  if(!fix.sigma){
    sig <- diag(D)
    sig[lower.tri(sig)] <- pv[((5*J)+1):length(pv)]
    sigm=Matrix::forceSymmetric(sig,uplo="L")
  }

  list(a=a,d=d,sigm=sigm)
}

P.Yi_2PL_3_rank <- function(a,#matrix 3 x J
                            d,#matrix 3 x J
                            theta, #matrix 3 x J  ability for each dimension on each item
                            yi# vector of J
){
  # yij is 1,...6
  response.pattern <- matrix(c(1,2,3,
                               1,3,2,
                               2,1,3,
                               2,3,1,
                               3,1,2,
                               3,2,1),ncol = 3,byrow = TRUE)
  z <- exp(a*theta-d) #3 x J
  sumz <- matrix(colSums(z)[col(z)],nrow = 3)
  red.sumz <- sumz-z
  p.a.abc <- z/sumz #three rows: za/sum(z),zb/sum(z),zc/sum(z); J columns
  p.b.bc <- z[response.pattern[,2],]/red.sumz[response.pattern[,1],]
  p.y <- p.a.abc[response.pattern[,1],]*p.b.bc
  p.y[p.y<1e-5] <- 1e-5
  p.y[matrix(c(yi,seq_len(length(yi))),ncol = 2)]
}

logLi_2PL_3_rank <- function(thetai, #vector of length D  ability for person i
                             a,#matrix 3 x J
                             d,#matrix 3 x J
                             yi,# vector of J
                             BID,
                             prior=TRUE,
                             sigma=NULL){
  ll <- sum(log(P.Yi_2PL_3_rank(a=a,d=d,theta=matrix(thetai[BID$Dim],nrow = 3),yi=yi)))
  lprior <- 0
  if(prior)
    lprior <- mvnfast::dmvn(X=matrix(thetai,nrow = 1),mu=rep(0,length(thetai)),sigma = sigma,log = TRUE)
  ll+lprior
}
logL_2PL_3_rank <- function(theta, #vector of length D  ability for person i
                            a,#matrix 3 x J
                            d,#matrix 3 x J
                            yi,# vector of J
                            BID,
                            sigma=NULL){
  # yij is 1,...6
  thetai <- matrix(theta[BID$Dim],nrow = 3)
  response.pattern <- matrix(c(1,2,3,
                               1,3,2,
                               2,1,3,
                               2,3,1,
                               3,1,2,
                               3,2,1),ncol = 3,byrow = TRUE)
  z <- exp(a*thetai-d) #3 x J
  sumz <- matrix(colSums(z)[col(z)],nrow = 3)
  red.sumz <- sumz-z
  p.a.abc <- z/sumz #three rows: za/sum(z),zb/sum(z),zc/sum(z); J columns
  p.b.bc <- z[response.pattern[,2],]/red.sumz[response.pattern[,1],]
  p.y <- p.a.abc[response.pattern[,1],]*p.b.bc[matrix(c(yi,seq_len(length(yi))),ncol = 2)]
  p.y[p.y<1e-5] <- 1e-5
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
