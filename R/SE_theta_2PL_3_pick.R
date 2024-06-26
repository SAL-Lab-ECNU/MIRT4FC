px_2PL_3_pick <- function(theta,a,d,xx=1){ #a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)
  z <- exp(a*theta-d)
  prod(z[xx])/sum(z)
}

grad_theta_2PL_3_pick <- function(theta,a,d,x=1){
  z <- exp(a*theta-d)
  ai <- a[1]
  aj <- a[2]
  ak <- a[3]
  # p.k.ijk <- p.k.ikj <- p.k.jki <- p.k.jik <- p.k.kij <- p.k.kji
  p.k <- z[3]/sum(z)
  # p.j.ijk <- p.j.ikj <- p.j.jki <- p.j.jik <- p.j.kij <- p.j.kji
  p.j <- z[2]/sum(z)
  # p.i.ijk <- p.i.ikj <- p.i.jki <- p.i.jik <- p.i.kij <- p.i.kji
  p.i <- z[1]/sum(z)
  if(identical(x,1)){
    g.i <- ai*p.i*(1-p.i)
    g.j <- -aj*p.i*p.j
    g.k <- -ak*p.i*p.k
  }else if(identical(x,2)){
    g.j <- aj*p.j*(1-p.j)
    g.k <- -ak*p.j*p.k
    g.i <- -ai*p.j*p.i
  }else if(identical(x,3)){
    g.k <- ak*p.k*(1-p.k)
    g.j <- -aj*p.k*p.j
    g.i <- -ai*p.k*p.i
  }

  c(g.i,g.j,g.k)
}

hes_theta_2PL_3_pick <- function(theta,a,d,x=1){
  h <- matrix(NA,3,3)
  z <- exp(a*theta-d)
  ai <- a[1]
  aj <- a[2]
  ak <- a[3]
  # p.k.ijk <- p.k.ikj <- p.k.jki <- p.k.jik <- p.k.kij <- p.k.kji
  p.k <- z[3]/sum(z)
  # p.j.ijk <- p.j.ikj <- p.j.jki <- p.j.jik <- p.j.kij <- p.j.kji
  p.j <- z[2]/sum(z)
  # p.i.ijk <- p.i.ikj <- p.i.jki <- p.i.jik <- p.i.kij <- p.i.kji
  p.i <- z[1]/sum(z)


  if(identical(x,1)){
    ii <- ai^2*p.i*(1-2*p.i)*(1-p.i)
    ij <- ji <- -aj*ai*p.i*p.j*(1-2*p.i)
    ik <- ki <- -ak*ai*p.i*p.k*(1-2*p.i)
    jj <- -aj^2*p.i*p.j*(1-2*p.j)
    kk <- -ak^2*p.i*p.k*(1-2*p.k)
    jk <- kj <- 2*aj*ak*p.i*p.j*p.k
  }else if(identical(x,2)){
    jj <- aj^2*p.j*(1-2*p.j)*(1-p.j)
    jk <- kj <- -ak*aj*p.j*p.k*(1-2*p.j)
    ji <- ij <- -ai*aj*p.j*p.i*(1-2*p.j)
    kk <- -ak^2*p.j*p.k*(1-2*p.k)
    ii <- -ai^2*p.j*p.i*(1-2*p.i)
    ki <- ik <- 2*ai*ak*p.i*p.j*p.k
  }else if(identical(x,3)){
    kk <- ak^2*p.k*(1-2*p.k)*(1-p.k)
    ki <- ik <- -ak*ai*p.i*p.k*(1-2*p.k)
    kj <- jk <- -ak*aj*p.j*p.k*(1-2*p.k)
    ii <- -ai^2*p.k*p.i*(1-2*p.i)
    jj <- -aj^2*p.k*p.j*(1-2*p.j)
    ij <- ji <- 2*ai*aj*p.i*p.j*p.k
  }
  h <- diag(c(ii,jj,kk))
  h[1,2] <- ij
  h[2,1] <- ji
  h[1,3] <- ik
  h[3,1] <- ki
  h[2,3] <- jk
  h[3,2] <- kj
  h
}


Ib_theta_2PL_3_pick <- function(theta,a,d){
  all.x <- list(1,2,3)
  ret <- 0
  for(x in 1:length(all.x)){
    g <- grad_theta_2PL_3_pick(theta,a,d,x=all.x[[x]])
    ret <- ret + tcrossprod(g,g)/px_2PL_3_pick(theta,a,d,xx=all.x[[x]])-hes_theta_2PL_3_pick(theta,a,d,x=all.x[[x]])
  }
  ret
}


test.info_theta_2PL_3_pick <- function(theta,item.par,BID){
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
    complete.info.b[dim.b,dim.b] <- info[[b]] <- Ib_theta_2PL_3_pick(theta.b,a.b,d.b)
    I <- I+complete.info.b
    comp.info[[b]] <- complete.info.b
  }
  se <- sqrt(diag(solve(Reduce("+",comp.info))))
  list(comp.info=comp.info,se =se, info.b=info,I=I)
}
