
logLi_2 <- function(thetai,
                  a,
                  d,
                  yi,
                  BID,
                  prior=TRUE,
                  sigma=NULL){
  ll <- sum(log(P.Yi_2(a=a,d=d,theta=matrix(thetai[BID$Dim],nrow = 2),yi=yi)))
  lprior <- 0
  if(prior)
    lprior <- mvnfast::dmvn(X=matrix(thetai,nrow = 1),mu=rep(0,length(thetai)),sigma = sigma,log = TRUE)
  ll+lprior
}


batch.var <- function(plist,n){
  phi.bar <- sapply(plist,rowMeans) # of par x n
  phi.hat <- rowMeans(Reduce(cbind,plist))
  rowMeans((phi.bar-phi.hat)^2)/(n-1)
}



Pj_2 <- function(j,theta,aj=NULL,dj=NULL,BID=NULL,item.par=NULL){
  response.pattern <- c(1,2)
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


  z <- exp(aj*t(thetaj)-dj) # 2 x N
  sumz <- matrix(colSums(z)[col(z)],nrow = 2)
  pj.a.ab <- z/sumz #three rows: za/sum(z),zb/sum(z); J columns
  pj.a.ab[response.pattern,,drop=FALSE]
}

#' data.sim_2
#'
#' @param item.par   a data.frame including a,d parameter
#' @param theta      matrix N x D
#' @param BID        block,item,dimension table
#'
#' @return N*J matrix.
#'
#'
#'
data.sim_2 <- function(item.par,theta,BID){

  Y <- matrix(NA,nrow = nrow(theta),ncol = max(BID$Block))
  for(j in unique(BID$Block)){
    pj <- Pj_2(j=j,theta=theta,BID=BID,item.par=item.par)

    Y[,j] <- t(apply(pj,2,function(p)sample(x=nrow(pj),size=1,prob = p)))

  }
  Y
}



kernel_2 <- function(Y,BID,a,d,sigm,theta,B=20,J,D,N,cor.matrix=FALSE,positive,fix.sigma=FALSE){

arm.sampler_2 <- function(thetai,p){
  logLi_2(thetai, a=a, d=d,yi=Y[i,], BID=BID,prior=TRUE,sigma=sigm)
}

objfun_2 <- function(x,thetaj,yj){
  -1*P.Yj_2(aj=(x[1:2]),dj=c(x[3],-x[3]),thetaj,yj)
}

#number of parameters
  npar <- ifelse(fix.sigma,3*J,3*J+D*(D-1)/2)
  #matrix for all parameters x iterations
  pmatrix <- matrix(NA,npar,B)
  # theta.list <- list()
  for(b in seq_len(B)){
    # sampling theta - Equation 4
      theta <- foreach (i=seq_len(N),.combine="rbind",.errorhandling = "stop",.inorder = FALSE,
                        .export = c("arm.sampler_2","logLi_2","P.Yi_2","a","d","Y","BID","sigm"))%dopar%{
                          armspp::arms_gibbs(n_samples = 1,previous = theta[i,],log_pdf = arm.sampler_2,
                                             lower = rep(-6,D),upper = rep(6,D),metropolis = fix.sigma)
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
      opt <- foreach (j=1:J,.combine="cbind",.errorhandling = "stop",.inorder = TRUE,
                      .export = c("objfun_2","P.Yj_2","a","d","Y","BID","theta"))%dopar%{
                       lo <- c(ifelse(positive[which(BID$Block==j)]>0,c(0,0),c(-Inf,-Inf)),-6)
                       up <- c(ifelse(positive[which(BID$Block==j)]>0,c(Inf,Inf),c(0,0)),6)
                       optim(par=c(a[,j],d[1,j]),fn = objfun_2,method = "L-BFGS-B",lower = lo, upper = up,
                              thetaj=theta[,BID$Dim[which(BID$Block==j)]],yj=Y[,j],control = list(maxit=20))$par
                      }



    a <- opt[1:2,]
    d <- rbind(opt[3,],-opt[3,])
    if(fix.sigma){
      pmatrix[,b] <- c(opt)
    }else{
      pmatrix[,b] <- c(c(opt),c(sigm[lower.tri(sigm)]))
    }


  }

  return(list(pmatrix=pmatrix,a=a,d=d,sigm=sigm,theta=theta))
}




P.Yj_2 <- function(aj,#vector of length 2
                 dj,#vector of length 2
                 thetaj,#matrix of N x 2
                 yj # vector of length N
){

  response.pattern <- c(1,2)
  z <- exp(aj*t(thetaj)-dj) # 2 x N
  sumz <- matrix(colSums(z)[col(z)],nrow = 2)
  p.a.ab <- z/sumz #three rows: za/sum(z),zb/sum(z),zc/sum(z); J columns
  p.y <- p.a.ab[response.pattern,]
  p <- p.y[matrix(c(yj,seq_len(length(yj))),ncol = 2)]
  p[p<1e-5] <- 1e-5
  p[p>1-1e-5] <- 1-1e-5
  sum(log(p))
}


P.Yi_2 <- function(a,#matrix 2 x J
                 d,#matrix 2 x J
                 theta, #matrix 2 x J  ability for each dimension on each item
                 yi# vector of J
){
  # yij is 1,2
  response.pattern <- c(1,2)
  z <- exp(a*theta-d) #3 x J
  sumz <- matrix(colSums(z)[col(z)],nrow = 2)
  p.a.ab <- z/sumz #three rows: za/sum(z),zb/sum(z); J columns
  p.y <- p.a.ab[response.pattern,]
  p.y[matrix(c(yi,seq_len(length(yi))),ncol = 2)]
}


param.vec.2.ads_2 <- function(pv,J,D,sigm=NULL,fix.sigma=FALSE){
  #pv: the vector of all parameters (a,d,sig)
  ad <- matrix(pv[1:(3*J)],nrow = 3)
  a <- ad[1:2,]
  d <- rbind(ad[3,],-ad[3,])
  if(!fix.sigma){
    sig <- diag(D)
    sig[lower.tri(sig)] <- pv[((3*J)+1):length(pv)]
    sigm=Matrix::forceSymmetric(sig,uplo="L")
  }

  list(a=a,d=d,sigm=sigm)
}
