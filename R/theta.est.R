#'theta.est
#'
#' @param Y          A # of subjects * # of blocks matrix; item responses.
#' @param BID        A # of statements * 3 matrix; item information,
#'     columns are "Block", "Item" and "Dimensions".
#' @param positive   A logical vector; indicating whether each statement is
#'     positive directional or not.
#' @param blocksize  A number; block size of FC(2/3/4).
#' @param res        A string; response format('pick'/'rank'/'mole'),
#'     pick-2(blocksize=2)/rank-2/mole-2 are equivalent, rank-3/mole-3 are equivalent.
#' @param M          A number; # of batch.
#' @param B          A number; # of iterations in each batch.
#' @param a          A vector; length = # of statements, initial alpha parameters.
#' @param d          A vector; length = # of statements, initial beta parameters.
#' @param item.par   A data frame; initial parameters for a and d.
#' @param sigma      A # of dimensions * # of dimensions matrix; initial sigma
#'     parameters.
#' @param theta      A # of subjects * # of dimensions matrix; initial theta parameters.
#' @param fix.sigma  Logical; TRUE if sigma is estimated.
#' @param burnin.maxitr   A number; max burn-in allowed.
#' @param maxitr     A number; max iterations allowed.
#' @param eps1       A number; stability criteria.
#' @param eps2       A number; convergence criterion.
#' @param frac1      A number; cutoffs for calculating Geweke z.
#' @param frac2      A number; cutoffs for calculating Geweke z.
#' @param cores      A number; number of cores.
#'
#' @return A list of data; including estimated a, d, and sigm parameters,
#'    total batch number, final chain size, burn-in size, time.
#'
#' @export theta.est
#' @examples
#' \donttest{
#' D <- 6
#' nitem.per.dim <- 10
#' nblock <- D * nitem.per.dim / 3
#' set.seed(123456)
#' item.par <- data.frame(a=seq_len(D*nitem.per.dim))
#' item.par <- within(item.par,{
#'   a <- runif(D*nitem.per.dim,0.7,3)
#'   b <- rnorm(D*nitem.per.dim)
#'   d <- a*b
#' })
#' BID <- data.frame(Block=rep(1:nblock,each=3),
#'                   Item=rep(1:3,nblock),
#'                   Dim=c(combn(D,3)[,sample(choose(D,3),nblock,replace = TRUE)]))
#' item.par$d <- c(t(aggregate(item.par$d,by=list(BID$Block),function(x)x-mean(x))[,-1]))
#' N <- 1000
#' v <- matrix(0.5,D,D)
#' diag(v) <- 1
#' theta <- mvnfast::rmvn(N,seq(-1,1,length.out = D),sigma = v)
#' Y <- data.sim(item.par,theta,BID,blocksize=3,res='rank')
#' x <- StEM(Y,BID,maxitr = 100,blocksize=3,res='rank',fix.sigma = TRUE,cores=1)
#' }
#' @importFrom stats cor optim rnorm
#' @importFrom utils install.packages installed.packages
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach registerDoSEQ %dopar%
#' @importFrom coda mcmc geweke.diag
#' @importFrom Matrix forceSymmetric
#' @importFrom parallel detectCores makeCluster stopCluster
theta.est <- function(Y,a=NULL,d=NULL,t=NULL,BID,sigma=NULL,prior=TRUE,blocksize=3,res='rank',model='2PL'){
  if(model=='2PL'){
    x <- theta.est_2PL(Y,a=a,d=d,BID,sigma=sigma,prior=prior,blocksize=blocksize,res=res,model=model)
  }else if(model=='TIRT'){
    x <- theta.est_TIRT(Y,a=a,d=d,BID,sigma=sigma,prior=prior,blocksize=blocksize,res=res,model=model)
  }else if(model=='ZG'){
    x <- theta.est_ZG(Y,a=a,t=t,BID,sigma=sigma,prior=prior,blocksize=blocksize,res=res,model=model)
  }else if(model=='GGUM'){
    x <- theta.est_GGUM(Y,a=a,d=d,t=t,BID,sigma=sigma,prior=prior,blocksize=blocksize,res=res,model=model)
  }
  return(x)
}

