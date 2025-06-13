#' Data.sim
#'
#'(item.par=item.par,theta=theta,BID=BID,blocksize=3,res='rank',model='2PL')
#' @param item.par   A data frame; parameters for a and d.
#' @param theta      A # of subjects * # of dimensions matrix; theta parameters.
#' @param BID        A # of statements * 3 matrix; item information,
#'     columns are "Block", "Item" and "Dimensions".
#' @param blocksize  A number; block size of FC(2/3/4).
#' @param res        A string; response format('pick'/'rank'/'mole'),
#'     pick-2(blocksize=2)/rank-2/mole-2 are equivalent, rank-3/mole-3 are equivalent.
#' @param model      A string; FC model('2PL'').
#'
#' @return      A # of subjects * # of block number matrix.
#' @export data.sim
#' @examples
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
#'

data.sim <- function(item.par=item.par,theta=theta,BID=BID,blocksize=3,res='rank',model='2PL'){
  if(model=='2PL'){
    Y <- data.sim_2PL(item.par=item.par,theta=theta,BID=BID,blocksize=blocksize,res=res)
  }
  return(Y)
}
