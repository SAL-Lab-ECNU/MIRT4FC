px_2PL_3_rank <- function(theta,a,d,x=c(1,2,3)){ #a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)
  z <- exp(a*theta-d)
  prod(z[x[-3]])/(sum(z)*sum(z[x[-1]]))
}

grad_theta_2PL_3_rank <- function(theta,a,d,x=c(1,2,3)){
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
  p.k.ki <- p.k.ik <- z[3]/sum(z[-2])
  p.i.ki <- p.i.ik <- z[1]/sum(z[-2])
  p.j.ji <- p.j.ij <- z[2]/sum(z[-3])
  p.i.ji <- p.i.ij <- z[1]/sum(z[-3])
  p.j.jk <- p.j.kj <- z[2]/sum(z[-1])
  p.k.jk <- p.k.kj <- z[3]/sum(z[-1])
  if(identical(x,c(1,2,3))){
    g.i <- ai*p.i*p.j
    g.j <- aj*p.i*p.j.jk*(p.k.jk-p.j)
    g.k <- -ak*p.i*p.j.jk*(p.k+p.k.jk)
  }else if(identical(x,c(2,3,1))){
    g.j <- aj*p.j*p.k
    g.k <- ak*p.j*p.k.ki*(p.i.ki-p.k)
    g.i <- -ai*p.j*p.k.ki*(p.i+p.i.ki)
  }else if(identical(x,c(1,3,2))){
    g.i <- ai*p.i*p.k
    g.k <- ak*p.i*p.k.kj*(p.j.kj-p.k)
    g.j <- -aj*p.i*p.k.kj*(p.j+p.j.kj)
  }else if(identical(x,c(2,1,3))){
    g.j <- aj*p.j*p.i
    g.i <- ai*p.j*p.i.ik*(p.k.ik-p.i)
    g.k <- -ak*p.j*p.i.ik*(p.k+p.k.ik)
  }else if(identical(x,c(3,1,2))){
    g.k <- ak*p.k*p.i
    g.i <- ai*p.k*p.i.ij*(p.j.ij-p.i)
    g.j <- -aj*p.k*p.i.ij*(p.j+p.j.ij)
  }else if(identical(x,c(3,2,1))){
    g.k <- ak*p.k*p.j
    g.j <- aj*p.k*p.j.ji*(p.i.ji-p.j)
    g.i <- -ai*p.k*p.j.ji*(p.i+p.i.ji)
  }

  c(g.i,g.j,g.k)
}

hes_theta_2PL_3_rank <- function(theta,a,d,x=c(1,2,3)){
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
  p.k.ki <- p.k.ik <- z[3]/sum(z[-2])
  p.i.ki <- p.i.ik <- z[1]/sum(z[-2])
  p.j.ji <- p.j.ij <- z[2]/sum(z[-3])
  p.i.ji <- p.i.ij <- z[1]/sum(z[-3])
  p.j.jk <- p.j.kj <- z[2]/sum(z[-1])
  p.k.jk <- p.k.kj <- z[3]/sum(z[-1])

  if(identical(x,c(1,2,3))){
    ii <- ai^2*p.i*p.j*(1-2*p.i)
    ij <- ji <- ai*aj*p.i*p.j*(1-2*p.j)
    ik <- ki <- -2*ai*ak*p.i*p.j*p.k
    jj <- -aj^2*p.i*p.j*p.j.jk*(1-2*p.j+p.k.jk)-aj^2*p.i*p.j.jk*p.k.jk*(p.j-p.k.jk+p.j.jk)
    kk <- ak^2*p.i*p.j.jk*((p.k+p.k.jk)^2+p.k^2-p.k-p.j.jk*p.k.jk)
    jk <- kj <- aj*ak*p.i*p.j.jk*(-1*p.k*p.k.jk+2*p.k*p.j-p.k.jk*p.k.jk+p.k.jk*p.j+p.j.jk*p.k.jk)
  }else if(identical(x,c(2,3,1))){
    jj <- aj^2*p.j*p.k*(1-2*p.j)
    jk <- kj <- aj*ak*p.j*p.k*(1-2*p.k)
    ji <- ij <- -2*aj*ai*p.j*p.k*p.i
    kk <- -ak^2*p.j*p.k*p.k.ki*(1-2*p.k+p.i.ki)-ak^2*p.j*p.k.ki*p.i.ki*(p.k-p.i.ki+p.k.ki)
    ii <- ai^2*p.j*p.k.ki*((p.i+p.i.ki)^2+p.i^2-p.i-p.k.ki*p.i.ki)
    ki <- ik <- ak*ai*p.j*p.k.ki*(-1*p.i*p.i.ki+2*p.i*p.k-p.i.ki*p.i.ki+p.i.ki*p.k+p.k.ki*p.i.ki)
  }else if(identical(x,c(1,3,2))){
    ii <- ai^2*p.i*p.k*(1-2*p.i)
    ik <- ki <- ai*ak*p.i*p.k*(1-2*p.k)
    ij <- ji <- -2*ai*aj*p.i*p.k*p.j
    kk <- -ak^2*p.i*p.k*p.k.kj*(1-2*p.k+p.j.kj)-ak^2*p.i*p.k.kj*p.j.kj*(p.k-p.j.kj+p.k.kj)
    jj <- aj^2*p.i*p.k.kj*((p.j+p.j.kj)^2+p.j^2-p.j-p.k.kj*p.j.kj)
    kj <- jk <- ak*aj*p.i*p.k.kj*(-1*p.j*p.j.kj+2*p.j*p.k-p.j.kj*p.j.kj+p.j.kj*p.k+p.k.kj*p.j.kj)
  }else if(identical(x,c(2,1,3))){
    jj <- aj^2*p.j*p.i*(1-2*p.j)
    ji <- ij <- aj*ai*p.j*p.i*(1-2*p.i)
    jk <- kj <- -2*aj*ak*p.j*p.i*p.k
    ii <- -ai^2*p.j*p.i*p.i.ik*(1-2*p.i+p.k.ik)-ai^2*p.j*p.i.ik*p.k.ik*(p.i-p.k.ik+p.i.ik)
    kk <- ak^2*p.j*p.i.ik*((p.k+p.k.ik)^2+p.k^2-p.k-p.i.ik*p.k.ik)
    ik <- ki <- ai*ak*p.j*p.i.ik*(-1*p.k*p.k.ik+2*p.k*p.i-p.k.ik*p.k.ik+p.k.ik*p.i+p.i.ik*p.k.ik)
  }else if(identical(x,c(3,1,2))){
    kk <- ak^2*p.k*p.i*(1-2*p.k)
    ki <- ik <- ak*ai*p.k*p.i*(1-2*p.i)
    kj <- jk <- -2*ak*aj*p.k*p.i*p.j
    ii <- -ai^2*p.k*p.i*p.i.ij*(1-2*p.i+p.j.ij)-ai^2*p.k*p.i.ij*p.j.ij*(p.i-p.j.ij+p.i.ij)
    jj <- aj^2*p.k*p.i.ij*((p.j+p.j.ij)^2+p.j^2-p.j-p.i.ij*p.j.ij)
    ij <- ji <- ai*aj*p.k*p.i.ij*(-1*p.j*p.j.ij+2*p.j*p.i-p.j.ij*p.j.ij+p.j.ij*p.i+p.i.ij*p.j.ij)
  }else if(identical(x,c(3,2,1))){
    kk <- ak^2*p.k*p.j*(1-2*p.k)
    kj <- jk <- ak*aj*p.k*p.j*(1-2*p.j)
    ki <- ik <- -2*ak*ai*p.k*p.j*p.i
    jj <- -aj^2*p.k*p.j*p.j.ji*(1-2*p.j+p.i.ji)-aj^2*p.k*p.j.ji*p.i.ji*(p.j-p.i.ji+p.j.ji)
    ii <- ai^2*p.k*p.j.ji*((p.i+p.i.ji)^2+p.i^2-p.i-p.j.ji*p.i.ji)
    ji <- ij <- aj*ai*p.k*p.j.ji*(-1*p.i*p.i.ji+2*p.i*p.j-p.i.ji*p.i.ji+p.i.ji*p.j+p.j.ji*p.i.ji)
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


Ib_theta_2PL_3_rank <- function(theta,a,d){
  all.x <- list(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
  ret <- 0
  for(x in 1:length(all.x)){
    g <- grad_theta_2PL_3_rank(theta,a,d,x=all.x[[x]])
    ret <- ret + tcrossprod(g,g)/px_2PL_3_rank(theta,a,d,x=all.x[[x]])-hes_theta_2PL_3_rank(theta,a,d,x=all.x[[x]])
  }
  ret
}


test.info_theta_2PL_3_rank <- function(theta,item.par,BID){
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
    complete.info.b[dim.b,dim.b] <- info[[b]] <- Ib_theta_2PL_3_rank(theta.b,a.b,d.b)
    I <- I+complete.info.b
    comp.info[[b]] <- complete.info.b
  }
  se <- sqrt(diag(solve(Reduce("+",comp.info))))
  list(comp.info=comp.info,se =se, info.b=info,I=I)
}
