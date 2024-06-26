px_2PL_4_pick <- function(theta,a,d,xx=1){ #a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)
  z <- exp(a*theta-d)
  prod(z[xx])/sum(z)
}

grad_theta_2PL_4_pick <- function(theta,a,d,x=1){
  z <- exp(a*theta-d)
  ai <- a[1]
  aj <- a[2]
  ak <- a[3]
  al <- a[4]
  p.l <- z[4]/sum(z)
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
    g.l <- -al*p.i*p.l
  }else if(identical(x,2)){
    g.j <- aj*p.j*(1-p.j)
    g.k <- -ak*p.j*p.k
    g.l <- -al*p.j*p.l
    g.i <- -ai*p.j*p.i
  }else if(identical(x,3)){
    g.k <- ak*p.k*(1-p.k)
    g.l <- -al*p.k*p.l
    g.i <- -ai*p.k*p.i
    g.j <- -aj*p.k*p.j
  }else if(identical(x,4)){
    g.l <- -al*p.l*(1-p.l)
    g.i <- -ai*p.l*p.i
    g.j <- -aj*p.l*p.j
    g.k <- -ak*p.l*p.k
  }

  c(g.i,g.j,g.k,g.l)
}

hes_theta_2PL_4_pick <- function(theta,a,d,x=1){
  h <- matrix(NA,4,4)
  z <- exp(a*theta-d)
  ai <- a[1]
  aj <- a[2]
  ak <- a[3]
  al <- a[4]
  p.l <- z[4]/sum(z)
  # p.k.ijk <- p.k.ikj <- p.k.jki <- p.k.jik <- p.k.kij <- p.k.kji
  p.k <- z[3]/sum(z)
  # p.j.ijk <- p.j.ikj <- p.j.jki <- p.j.jik <- p.j.kij <- p.j.kji
  p.j <- z[2]/sum(z)
  # p.i.ijk <- p.i.ikj <- p.i.jki <- p.i.jik <- p.i.kij <- p.i.kji
  p.i <- z[1]/sum(z)


  if(identical(x,1)){
    ii <- ai^2*p.i*(1-2*p.i)
    ij <- ji <- -aj*ai*p.i*p.j*(1-2*p.i)
    ik <- ki <- -ak*ai*p.i*p.k*(1-2*p.i)
    il <- li <- -al*ai*p.i*p.l*(1-2*p.i)
    jj <- -aj^2*p.i*p.j*(1-2*p.j)
    jk <- kj <- 2*aj*ak*p.i*p.j*p.k
    jl <- lj <- 2*aj*al*p.i*p.j*p.l
    kk <- -ak^2*p.i*p.k*(1-2*p.k)
    kl <-  lk <- 2*al*ak*p.i*p.l*p.k
    ll <- -al^2*p.i*p.l*(1-2*p.l)
  }else if(identical(x,2)){
    jj <- aj^2*p.j*(1-2*p.j)
    jk <- kj <- -ak*aj*p.j*p.k*(1-2*p.j)
    jl <- lj <- -al*aj*p.j*p.l*(1-2*p.j)
    ji <- ij <- -ai*aj*p.j*p.i*(1-2*p.j)
    kk <- -ak^2*p.j*p.k*(1-2*p.k)
    kl <- lk <- 2*al*ak*p.l*p.j*p.k
    ki <- ik <- 2*ai*ak*p.i*p.j*p.k
    ll <- -al^2*p.i*p.l*(1-2*p.l)
    li <- il <- 2*ai*al*p.i*p.j*p.l
    ii <- -ai^2*p.j*p.i*(1-2*p.i)
  }else if(identical(x,3)){
  kk <- ak^2*p.k*(1-2*p.k)
  kl <- lk <- -al*ak*p.k*p.l*(1-2*p.k)
  ki <- ik <- -ai*ak*p.k*p.i*(1-2*p.k)
  kj <- jk <- -aj*ak*p.k*p.j*(1-2*p.k)
  ll <- -al^2*p.k*p.l*(1-2*p.l)
  li <- il <- 2*ai*al*p.i*p.k*p.l
  lj <- jl <- 2*aj*al*p.j*p.k*p.l
  ii <- -ai^2*p.j*p.i*(1-2*p.i)
  ij <- ji <- 2*aj*ai*p.j*p.k*p.i
  jj <- -aj^2*p.k*p.j*(1-2*p.j)
}else if(identical(x,4)){
  ll <- al^2*p.l*(1-2*p.l)
  li <- il <- -ai*al*p.l*p.i*(1-2*p.l)
  lj <- jl <- -aj*al*p.l*p.j*(1-2*p.l)
  lk <- kl <- -ak*al*p.l*p.k*(1-2*p.l)
  ii <- -ai^2*p.l*p.i*(1-2*p.i)
  ij <- ji <- 2*aj*ai*p.j*p.l*p.i
  ik <- ki <- 2*ak*ai*p.k*p.l*p.i
  jj <- -aj^2*p.k*p.j*(1-2*p.j)
  jk <- kj <- 2*ak*aj*p.k*p.l*p.j
  kk <- -ak^2*p.l*p.k*(1-2*p.k)
}
  h <- diag(c(ii,jj,kk,ll))
  h[1,2] <- ij
  h[2,1] <- ji
  h[1,3] <- ik
  h[3,1] <- ki
  h[1,4] <- il
  h[4,1] <- li
  h[2,3] <- jk
  h[3,2] <- kj
  h[2,4] <- jl
  h[4,2] <- lj
  h[3,4] <- kl
  h[4,3] <- lk
  h
}


Ib_theta_2PL_4_pick <- function(theta,a,d){
  all.x <- list(1,2,3,4)
  ret <- 0
  for(x in 1:length(all.x)){
    g <- grad_theta_2PL_4_pick(theta,a,d,x=all.x[[x]])
    ret <- ret + tcrossprod(g,g)/px_2PL_4_pick(theta,a,d,xx=all.x[[x]])-hes_theta_2PL_4_pick(theta,a,d,x=all.x[[x]])
  }
  ret
}


test.info_theta_2PL_4_pick <- function(theta,item.par,BID){
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
    complete.info.b[dim.b,dim.b] <- info[[b]] <- Ib_theta_2PL_4_pick(theta.b,a.b,d.b)
    I <- I+complete.info.b
    comp.info[[b]] <- complete.info.b
  }
  se <- sqrt(diag(solve(Reduce("+",comp.info))))
  list(comp.info=comp.info,se =se, info.b=info,I=I)
}
