px_TIRT_3_pick <- function(theta,a,d,xx=1){ #a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(ai,aj,ak)
  alphatheta=a*theta
  p12=pnorm(alphatheta[1]-alphatheta[2]-d[1]+d[2])
  p13=pnorm(alphatheta[1]-alphatheta[3]-d[1]+d[3])
  p21=1-p12
  p23=pnorm(alphatheta[2]-alphatheta[3]-d[2]+d[3])
  p31=1-p13
  p32=1-p23
  pa=p12*p13
  pb=p21*p23
  pc=p31*p32
  pall=pa+pb+pc
  if(identical(xx,1)){       #顺序为1,2,3时的一阶导，下方同理
    L=pa/pall
  }else if(identical(xx,2)){
    L=pb/pall
  }else if(identical(xx,3)){
    L=pc/pall
  }
  L
}

grad_theta_TIRT_3_pick <- function(theta,a,d,x=1){
  ai <- a[1]
  aj <- a[2]
  ak <- a[3]
  alphatheta=a*theta
  p12=pnorm(alphatheta[1]-alphatheta[2]-d[1]+d[2])
  p13=pnorm(alphatheta[1]-alphatheta[3]-d[1]+d[3])
  p21=1-p12
  p23=pnorm(alphatheta[2]-alphatheta[3]-d[2]+d[3])
  p31=1-p13
  p32=1-p23
  pa=p12*p13
  pb=p21*p23
  pc=p31*p32
  pall=pa+pb+pc
  if(identical(x,1)){
    g.i <- (ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13+p12*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))/pall-
      pa/pall^2*((ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13+p12*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+(-ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23)+(-ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32))
    g.j <- (-aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13)/pall-
      pa/pall^2*((-aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13)+(aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p23+p21*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+(-aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)*p31))
    g.k <- (-p12*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))/pall-
      pa/pall^2*((-p12*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+(-p21*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+(p32*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)+p31*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)))

  }else if(identical(x,2)){
    g.i <- (-ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23)/pall-
      pb/pall^2*((ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13+p12*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+(-ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23)+(-ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32))
    g.j <- (aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p23+p21*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))/pall-
      pb/pall^2*((-aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13)+(aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p23+p21*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+(-aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)*p31))
    g.k <- (-p21*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))/pall-
      pb/pall^2*((-p12*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+(-p21*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+(p32*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)+p31*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)))

  }else if(identical(x,3)){
    g.i <- (-ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32)/pall-
      pc/pall^2*((ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13+p12*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+(-ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23)+(-ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32))
    g.j <- (-aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)*p31)/pall-
      pc/pall^2*((-aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13)+(aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p23+p21*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+(-aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)*p31))
    g.k <- (p32*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)+p31*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))/pall-
      pc/pall^2*((-p12*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+(-p21*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+(p32*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)+p31*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)))

  }

  c(g.i,g.j,g.k)
}



Ib_theta_TIRT_3_pick <- function(theta,a,d){
  all.x <- list(1,2,3)
  ret <- 0
  for(x in 1:length(all.x)){
    g <- grad_theta_TIRT_3_pick(theta,a,d,x=all.x[[x]])
    ret <- ret + tcrossprod(g,g)/px_TIRT_3_pick(theta,a,d,xx=all.x[[x]])
  }
  ret
}


test.info_theta_TIRT_3_pick <- function(theta,item.par,BID){
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
    complete.info.b[dim.b,dim.b] <- info[[b]] <- Ib_theta_TIRT_3_pick(theta.b,a.b,d.b)
    I <- I+complete.info.b
    comp.info[[b]] <- complete.info.b
  }
  se <- sqrt(diag(solve(Reduce("+",comp.info))))
  list(comp.info=comp.info,se =se, info.b=info,I=I)
}
