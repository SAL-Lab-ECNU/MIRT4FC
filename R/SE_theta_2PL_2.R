px_2PL_2 <- function(theta,a,d,xx=1){ #a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)
  z <- exp(a*theta-d)
  prod(z[xx])/sum(z)
}

grad_theta_2PL_2 <- function(theta,a,d,x=1){
  z <- exp(a*theta-d)
  ai <- a[1]
  aj <- a[2]


  # p.j.ijk <- p.j.ikj <- p.j.jki <- p.j.jik <- p.j.kij <- p.j.kji
  p.j <- z[2]/sum(z)
  # p.i.ijk <- p.i.ikj <- p.i.jki <- p.i.jik <- p.i.kij <- p.i.kji
  p.i <- z[1]/sum(z)
  if(identical(x,1)){
    g.i <- ai*p.i*p.j
    g.j <- -aj*p.i*p.j
  }else if(identical(x,2)){
    g.j <- aj*p.i*p.j
    g.i <- -ai*p.j*p.i
  }

  c(g.i,g.j)
}

hes_theta_2PL_2 <- function(theta,a,d,x=1){
  h <- matrix(NA,2,2)
  z <- exp(a*theta-d)
  ai <- a[1]
  aj <- a[2]
  # p.j.ijk <- p.j.ikj <- p.j.jki <- p.j.jik <- p.j.kij <- p.j.kji
  p.j <- z[2]/sum(z)
  # p.i.ijk <- p.i.ikj <- p.i.jki <- p.i.jik <- p.i.kij <- p.i.kji
  p.i <- z[1]/sum(z)


  if(identical(x,1)){
    ii <- ai^2*p.i*p.j*(p.j-p.i)
    ij <- ji <- -aj*ai*p.i*p.j*(p.j-p.i)
    jj <- aj^2*p.i*p.j*(p.j-p.i)
  }else if(identical(x,2)){
    jj <- aj^2*p.i*p.j*(p.i-p.j)
    ji <- ij <- -aj*ai*p.i*p.j*(p.i-p.j)
    ii <- ai^2*p.i*p.j*(p.i-p.j)
  }
  h <- diag(c(ii,jj))
  h[1,2] <- ij
  h[2,1] <- ji
  h
}


Ib_theta_2PL_2 <- function(theta,a,d){
  all.x <- list(1,2)
  ret <- 0
  for(x in 1:length(all.x)){
    g <- grad_theta_2PL_2(theta,a,d,x=all.x[[x]])
    ret <- ret + tcrossprod(g,g)/px_2PL_2(theta,a,d,xx=all.x[[x]])-hes_theta_2PL_2(theta,a,d,x=all.x[[x]])
  }
  ret
}


test.info_theta_2PL_2 <- function(theta,item.par,BID){
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
    complete.info.b[dim.b,dim.b] <- info[[b]] <- Ib_theta_2PL_2(theta.b,a.b,d.b)
    I <- I+complete.info.b
    comp.info[[b]] <- complete.info.b
  }
  se <- sqrt(diag(solve(Reduce("+",comp.info))))
  list(comp.info=comp.info,se =se, info.b=info,I=I)
}
