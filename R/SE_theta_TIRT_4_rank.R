px_TIRT_4_rank <- function(theta,a,d,xx=c(1,2,3,4)){ #a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)
  alphatheta=a*theta
  p12=pnorm(alphatheta[1]-alphatheta[2]-d[1]+d[2])
  p13=pnorm(alphatheta[1]-alphatheta[3]-d[1]+d[3])
  p14=pnorm(alphatheta[1]-alphatheta[4]-d[1]+d[4])
  p21=1-p12
  p23=pnorm(alphatheta[2]-alphatheta[3]-d[2]+d[3])
  p24=pnorm(alphatheta[2]-alphatheta[4]-d[2]+d[4])
  p31=1-p13
  p32=1-p23
  p34=pnorm(alphatheta[3]-alphatheta[4]-d[3]+d[4])
  p41=1-p14
  p42=1-p24
  p43=1-p34
  pabcd=p12*p13*p14*p23*p24*p34
  pabdc=p12*p13*p14*p23*p24*p43
  pacbd=p12*p13*p14*p32*p24*p34
  pacdb=p12*p13*p14*p32*p42*p34
  padbc=p12*p13*p14*p23*p42*p43
  padcb=p12*p13*p14*p32*p42*p43

  pbacd=p21*p13*p14*p23*p24*p34
  pbadc=p21*p13*p14*p23*p24*p43
  pbcad=p21*p31*p14*p23*p24*p34
  pbcda=p21*p31*p41*p23*p24*p34
  pbdac=p21*p13*p41*p23*p24*p43
  pbdca=p21*p31*p41*p23*p24*p43

  pcabd=p12*p31*p14*p32*p24*p34
  pcadb=p12*p31*p14*p32*p42*p34
  pcbad=p21*p31*p14*p32*p24*p34
  pcbda=p21*p31*p41*p32*p24*p34
  pcdab=p12*p31*p41*p32*p42*p34
  pcdba=p21*p31*p41*p32*p42*p34

  pdabc=p12*p13*p41*p23*p42*p43
  pdacb=p12*p13*p41*p32*p42*p43
  pdbac=p21*p13*p41*p23*p42*p43
  pdbca=p21*p31*p41*p23*p42*p43
  pdcab=p12*p31*p41*p32*p42*p43
  pdcba=p21*p31*p41*p32*p42*p43
  pall=pabcd+pabdc+pacbd+pacdb+padbc+padcb+pbacd+pbadc+pbcad+pbcda+pbdac+pbdca+
    pcabd+pcadb+pcbad+pcbda+pcdab+pcdba+pdabc+pdacb+pdbac+pdbca+pdcab+pdcba
  if(identical(xx,c(1,2,3,4))){       #顺序为1,2,3时的一阶导，下方同理
    L=pabcd/pall
  }else if(identical(xx,c(1,2,4,3))){
    L=pabdc/pall
  }else if(identical(xx,c(1,3,2,4))){
    L=pacbd/pall
  }else if(identical(xx,c(1,3,4,2))){
    L=pacdb/pall
  }else if(identical(xx,c(1,4,2,3))){
    L=padbc/pall
  }else if(identical(xx,c(1,4,3,2))){
    L=padcb/pall
  }else if(identical(xx,c(2,1,3,4))){
    L=pbacd/pall
  }else if(identical(xx,c(2,1,4,3))){
    L=pbadc/pall
  }else if(identical(xx,c(2,3,1,4))){
    L=pbcad/pall
  }else if(identical(xx,c(2,3,4,1))){
    L=pbcda/pall
  }else if(identical(xx,c(2,4,1,3))){
    L=pbdac/pall
  }else if(identical(xx,c(2,4,3,1))){
    L=pbdca/pall
  }else if(identical(xx,c(3,1,2,4))){
    L=pcabd/pall
  }else if(identical(xx,c(3,1,4,2))){
    L=pcadb/pall
  }else if(identical(xx,c(3,2,1,4))){
    L=pcbad/pall
  }else if(identical(xx,c(3,2,4,1))){
    L=pcbda/pall
  }else if(identical(xx,c(3,4,1,2))){
    L=pcdab/pall
  }else if(identical(xx,c(3,4,2,1))){
    L=pcdba/pall
  }else if(identical(xx,c(4,1,2,3))){
    L=pdabc/pall
  }else if(identical(xx,c(4,1,3,2))){
    L=pdacb/pall
  }else if(identical(xx,c(4,2,1,3))){
    L=pdbac/pall
  }else if(identical(xx,c(4,2,3,1))){
    L=pdbca/pall
  }else if(identical(xx,c(4,3,1,2))){
    L=pdcab/pall
  }else if(identical(xx,c(4,3,2,1))){
    L=pdcba/pall
  }


  L
}

grad_theta_TIRT_4_rank <- function(theta,a,d,x=c(1,2,3,4),h=1e-5){
  g.i <- (px_TIRT_4_rank(theta+c(h,0,0,0),a,d,x)-px_TIRT_4_rank(theta,a,d,x))/h
  g.j <- (px_TIRT_4_rank(theta+c(0,h,0,0),a,d,x)-px_TIRT_4_rank(theta,a,d,x))/h
  g.k <- (px_TIRT_4_rank(theta+c(0,0,h,0),a,d,x)-px_TIRT_4_rank(theta,a,d,x))/h
  g.l <- (px_TIRT_4_rank(theta+c(0,0,0,h),a,d,x)-px_TIRT_4_rank(theta,a,d,x))/h
  c(g.i,g.j,g.k,g.l)
}



Ib_theta_TIRT_4_rank <- function(theta,a,d){
  all.x <- list(c(1,2,3,4),c(1,2,4,3),c(1,3,2,4),c(1,3,4,2),c(1,4,2,3),c(1,4,3,2),
                c(2,1,3,4),c(2,1,4,3),c(2,3,1,4),c(2,3,4,1),c(2,4,1,3),c(2,4,3,1),
                c(3,1,2,4),c(3,1,4,2),c(3,2,1,4),c(3,2,4,1),c(3,4,1,2),c(3,4,2,1),
                c(4,1,2,3),c(4,1,3,2),c(4,2,1,3),c(4,2,3,1),c(4,3,1,2),c(4,3,2,1))
  ret <- 0
  for(x in 1:length(all.x)){
    g <- grad_theta_TIRT_4_rank(theta,a,d,x=all.x[[x]])
    ret <- ret + tcrossprod(g,g)/px_TIRT_4_rank(theta,a,d,xx=all.x[[x]])
  }
  ret
}


test.info_theta_TIRT_4_rank <- function(theta,item.par,BID){
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
    complete.info.b[dim.b,dim.b] <- info[[b]] <- Ib_theta_TIRT_4_rank(theta.b,a.b,d.b)
    I <- I+complete.info.b
    comp.info[[b]] <- complete.info.b
  }
  se <- sqrt(diag(solve(Reduce("+",comp.info))))
  list(comp.info=comp.info,se =se, info.b=info,I=I)
}
