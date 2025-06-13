#'theta.est
#'(Y,a=NULL,d=NULL,BID,sigma=NULL,prior=TRUE,blocksize=3,res='rank',model='2PL')
#' @param Y          A # of subjects * # of blocks matrix; item responses.
#' @param BID        A # of statements * 3 matrix; item information,
#'     columns are "Block", "Item" and "Dimensions".
#' @param blocksize  A number; block size of FC(2/3/4).
#' @param res        A string; response format('pick'/'rank'/'mole'),
#'     pick-2(blocksize=2)/rank-2/mole-2 are equivalent, rank-3/mole-3 are equivalent.
#' @param a          A vector; length = # of statements, initial alpha parameters.
#' @param d          A vector; length = # of statements, initial beta parameters.
#' @param sigma      A # of dimensions * # of dimensions matrix; initial sigma
#'     parameters.
#' @param prior      Logical; TRUE if prior is added
#' @param model      A string; FC model('2PL'/'TIRT').
#'
#' @return A list of data; including estimated a, d, and sigm parameters,
#'    total batch number, final chain size, burn-in size, time.
#' @examples
#' D <- 6
#'nitem.per.dim <- 10
#'nblock <- D * nitem.per.dim / 3
#'set.seed(123456)
#'
#'item.par <- data.frame(a=seq_len(D*nitem.per.dim))
#'item.par <- within(item.par,{
#'  a <- runif(D*nitem.per.dim,0.7,3)
#'  b <- rnorm(D*nitem.per.dim)
#'  d <- a*b
#'})
#'a=matrix(item.par$a,nrow=3,ncol = 20)
#'d=matrix(item.par$d,nrow=3,ncol = 20)
#'
#'BID <- data.frame(Block=rep(1:nblock,each=3),
#'                  Item=rep(1:3,nblock),
#'                  Dim=c(combn(D,3)[,sample(choose(D,3),nblock,replace = TRUE)]))
#'item.par$d <- c(t(aggregate(item.par$d,by=list(BID$Block),function(x)x-mean(x))[,-1]))
#'N <- 1000
#'
#'v <- matrix(0,D,D)
#'diag(v) <- 1
#'eigen(v)$values
#'theta <- mvnfast::rmvn(N,rep(0,each = D),sigma = v)
#'Y <- data.sim(item.par,theta,BID,blocksize = 3,res = "rank")
#'
#'thetaest <- theta.est(Y,a,d,BID=BID,sigma=v,prior=TRUE,blocksize=3,res='rank',model='2PL')
#'
#' @export theta.est
#' @importFrom stats cor optim rnorm
#' @importFrom utils install.packages installed.packages
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach registerDoSEQ %dopar%
#' @importFrom coda mcmc geweke.diag
#' @importFrom Matrix forceSymmetric
#' @importFrom parallel detectCores makeCluster stopCluster
theta.est <- function(Y,a=NULL,d=NULL,BID,sigma=NULL,prior=TRUE,blocksize=3,res='rank',model='2PL'){
  if(model=='2PL'){
    x <- theta.est_2PL(Y,a=a,d=d,BID,sigma=sigma,prior=prior,blocksize=blocksize,res=res,model=model)
  }
  return(x)
}

