#' logLi_3_rank
#'
#' @param thetai  A vector; D(dimension) ability for person i.
#' @param a       A matrix 3 x J.
#' @param d       A matrix 3 x J.
#' @param yi      A vector; J(block number) response for person i.
#' @param BID     A # of statements x 3 matrix; item information,
#'     columns are "Block", "Item" and "Dimensions".
#' @param prior   Logical.
#' @param sigma   Matrix D x D or NULL.
#'
#' @return a number.
#'
#'
#'
logLi_3_rank <- function(thetai, #vector of length D  ability for person i
                    a,#matrix 3 x J
                    d,#matrix 3 x J
                    yi,# vector of J
                    BID,
                    prior=TRUE,
                    sigma=NULL){
  ll <- sum(log(P.Yi_3_rank(a=a,d=d,theta=matrix(thetai[BID$Dim],nrow = 3),yi=yi)))
  lprior <- 0
  if(prior)
    lprior <- mvnfast::dmvn(X=matrix(thetai,nrow = 1),mu=rep(0,length(thetai)),sigma = sigma,log = TRUE)
  ll+lprior
}

#' Pj_3_rank
#' @description Calculate probabilities of choosing each option out of
#'              6 alternatives of block j by all N persons: 6 x N matrix
#' @param j       block number
#' @param theta   matrix N x D
#' @param aj      matrix 1 x 3
#' @param dj      matrix 1 x 3
#' @param BID     block,item,dimension table
#' @param item.par  a data.frame including a,d parameter
#'
#' @return a 3*N matrix.
#'
#'
#'
Pj_3_rank <- function(j,theta,aj=NULL,dj=NULL,BID=NULL,item.par=NULL){
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
#' data.sim_3_rank
#'
#' @param item.par   a data.frame including a,d parameter
#' @param theta      matrix N x D
#' @param BID        block,item,dimension table
#'
#' @return N*J matrix.
#'
#'
#'
data.sim_3_rank <- function(item.par,theta,BID){

  Y <- matrix(NA,nrow = nrow(theta),ncol = max(BID$Block))
  for(j in unique(BID$Block)){
    pj <- Pj_3_rank(j=j,theta=theta,BID=BID,item.par=item.par)

    Y[,j] <- t(apply(pj,2,function(p)sample(x=nrow(pj),size=1,prob = p)))

  }
  Y
}

#' kernel_3_rank
#'
#' @param Y          a # of subjects x # of blocks matrix of item responses
#' @param BID        a # of statements x 3 matrix of item information; Columns are "Block", "Item" and "Dimensions"
#' @param a          initial alpha parameters ; a vector with length = # of statements
#' @param d          initial beta parameters ; a vector with length = # of statements
#' @param sigm       initial sigma parameters; a matrix with # of dimensions x # of dimensions
#' @param theta      initial theta parameters; a matrix with # of subjects x # of dimensions
#' @param B          of iterations in each batch
#' @param J          number of blocks
#' @param D          number of dimensions
#' @param N          number of participants
#' @param cor.matrix logical; use correlation matrix as the estimate if cor.matrix = TRUE
#' @param positive   a logical vector indicating whether each statement is positive directional or not
#' @param fix.sigma  logical; TRUE if sigma is not estimated
#'
#' @return A list of data.
#'
#'
#'
kernel_3_rank <- function(Y,BID,a,d,sigm,theta,B=20,J,D,N,cor.matrix=FALSE,positive,fix.sigma=FALSE){
# posterior for sampling theta - see Equation 4
arm.sampler_3_rank <- function(thetai,p){
  logLi_3_rank(thetai, a=a, d=d,yi=Y[i,], BID=BID,prior=TRUE,sigma=sigm)
}
# function to be optimized for estimating item parameters - see Equation 8
objfun_3_rank <- function(x,thetaj,yj){
  -1*P.Yj_3_rank(aj=(x[1:3]),dj=c(x[4:5],-sum(x[4:5])),thetaj,yj)
}

  #number of parameters
  npar <- ifelse(fix.sigma,5*J,5*J+D*(D-1)/2)
  #matrix for all parameters x iterations
  pmatrix <- matrix(NA,npar,B)
  # theta.list <- list()
  for(b in seq_len(B)){
    # sampling theta - Equation 4
    theta <- foreach (i=seq_len(N),.combine="rbind",.errorhandling = "stop",.inorder = TRUE,
                      .export = c("arm.sampler_3_rank","logLi_3_rank","P.Yi_3_rank","a","d","Y","BID","sigm"))%dopar%{

                        armspp::arms_gibbs(n_samples = 1,previous = theta[i,],log_pdf = arm.sampler_3_rank,
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
                    .export = c("objfun_3_rank","P.Yj_3_rank","a","d","Y","BID","theta"))%dopar%{
                      lo <- c(ifelse(positive[which(BID$Block==j)]>0,c(0,0,0),c(-Inf,-Inf,-Inf)),-6,-6)
                      up <- c(ifelse(positive[which(BID$Block==j)]>0,c(Inf,Inf,Inf),c(0,0,0)),6,6)
                      optim(par=c(a[,j],d[1:2,j]),fn = objfun_3_rank,method = "L-BFGS-B",lower = lo, upper = up,
                            thetaj=theta[,BID$Dim[which(BID$Block==j)]],yj=Y[,j],control = list(maxit=20))$par
                    }



    a <- opt[1:3,]
    d <- rbind(opt[4:5,],-1*colSums(opt[4:5,]))
    if(fix.sigma){
      pmatrix[,b] <- c(opt)
    }else{
      pmatrix[,b] <- c(c(opt),c(sigm[lower.tri(sigm)]))
    }


  }

  return(list(pmatrix=pmatrix,a=a,d=d,sigm=sigm,theta=theta))
}


#' P.Yj_3_rank
#' @description sum_i log P(yij|theta_i,aj,dj)
#'
#' @param aj  vector of length 3
#' @param dj  vector of length 3
#' @param thetaj matrix of N x 3
#' @param yj  vector of length N
#'
#' @return a number.
#'
#'
#'
P.Yj_3_rank <- function(aj,#vector of length 3
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
#' P.Yi_3_rank
#' @description probability of yi - length of J p(yi|theta,...) ---- verified
#'
#' @param a       matrix 3 x J
#' @param d       matrix 3 x J
#' @param theta   matrix 3 x J  ability for each dimension on each item
#' @param yi      vector of J
#'
#' @return A 3*J matrix.
#'
#'
#'
P.Yi_3_rank <- function(a,#matrix 3 x J
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
#' param.vec.2.ads_3_rank
#' @description function used to change the format of the results
#'
#' @param pv   the vector of all parameters (a,d,sig)
#' @param J    number of blocks
#' @param D    number of dimension
#' @param sigm       initial sigma parameters; a matrix with # of dimensions x # of dimensions
#' @param fix.sigma  logical; TRUE if sigma is not estimated
#'
#' @return A list of data.
#'
#'
#'
param.vec.2.ads_3_rank <- function(pv,J,D,sigm=NULL,fix.sigma=FALSE){
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
