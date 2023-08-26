#' Data.sim
#'
#' @param item.par   A data frame; parameters for a and d.
#' @param theta      A # of subjects x # of dimensions matrix; theta parameters.
#' @param BID        A # of statements x 3 matrix; item information,
#'     columns are "Block", "Item" and "Dimensions".
#' @param blocksize  A number; block size of FC(2/3/4).
#' @param res        A string; response format('pick'/'rank'/'mole').
#'
#' @return      A # of subjects x # of block number matrix.
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

data.sim <- function(item.par=item.par,theta=theta,BID=BID,blocksize=3,res='rank'){
  if(blocksize==2){
    Y <- data.sim_2(item.par=item.par,theta=theta,BID=BID)
  }else if(blocksize==3){
    if(res=='pick'){
      Y <- data.sim_3_pick(item.par=item.par,theta=theta,BID=BID)
    }else if(res=='rank'|res=='mole'){
      Y <- data.sim_3_rank(item.par=item.par,theta=theta,BID=BID)
    }
  }else if(blocksize==4){
    if(res=='pick'){
      Y <- data.sim_4_pick(item.par=item.par,theta=theta,BID=BID)
    }else if(res=='rank'){
      Y <- data.sim_4_rank(item.par=item.par,theta=theta,BID=BID)
    }else if(res=='mole'){
      Y <- data.sim_4_mole(item.par=item.par,theta=theta,BID=BID)
    }
  }
}
