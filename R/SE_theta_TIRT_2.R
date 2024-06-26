px_TIRT_2 <- function(theta,a,d,xx=1){ #a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(ai,aj,thetak)
  alphatheta=a*theta
  p12=pnorm(alphatheta[1]-alphatheta[2]-d[1]+d[2])
  p21=1-p12
  if(identical(xx,1)){       #顺序为1,2,3时的一阶导，下方同理
    L=p12
  }else if(identical(xx,2)){
    L=p21
  }
  L
}

grad_theta_TIRT_2 <- function(theta,a,d,x=1){
  alphatheta=a*theta
  ai <- a[1]
  aj <- a[2]


  if(identical(x,1)){
    g.i <- ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)
    g.j <- -aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)
  }else if(identical(x,2)){
    g.i <- -ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)
    g.j <- aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)
  }

  c(g.i,g.j)
}



Ib_theta_TIRT_2 <- function(theta,a,d){
  all.x <- list(1,2)
  ret <- 0
  for(x in 1:length(all.x)){
    g <- grad_theta_TIRT_2(theta,a,d,x=all.x[[x]])
    ret <- ret + tcrossprod(g,g)/px_TIRT_2(theta,a,d,xx=all.x[[x]])
  }
  ret
}


test.info_theta_TIRT_2 <- function(theta,item.par,BID){
  D <- max(BID$Dim) #number of dims
  comp.info <- info <- list()
  I <- diag(rep(0,D))
  for(b in unique(BID$Block)){
    complete.info.b <- diag(rep(0,D))
    items <- which(BID$Block==b)
    dim.b <- BID$Dim[items]
    a.b <- item.par$a[items]
    d.b <- item.par$d[items]
    theta.b <- theta[dim.b]
    complete.info.b[dim.b,dim.b] <- info[[b]] <- Ib_theta_TIRT_2(theta.b,a.b,d.b)
    I <- I+complete.info.b
    comp.info[[b]] <- complete.info.b
  }
  se <- sqrt(diag(solve(Reduce("+",comp.info))))
  list(comp.info=comp.info,se =se, info.b=info,I=I)
}
