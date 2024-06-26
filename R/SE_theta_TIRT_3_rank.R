px_TIRT_3_rank <- function(theta,a,d,xx=c(1,2,3)){ #a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(ai,aj,ak)
  alphatheta=a*theta
  p12=pnorm(alphatheta[1]-alphatheta[2]-d[1]+d[2])
  p13=pnorm(alphatheta[1]-alphatheta[3]-d[1]+d[3])
  p21=1-p12
  p23=pnorm(alphatheta[2]-alphatheta[3]-d[2]+d[3])
  p31=1-p13
  p32=1-p23
  pabc=p12*p13*p23
  pacb=p12*p13*p32
  pbac=p21*p23*p13
  pbca=p21*p23*p31
  pcab=p31*p32*p12
  pcba=p31*p32*p21
  pall=pabc+pacb+pbac+pbca+pcab+pcba
  if(identical(xx,c(1,2,3))){       #顺序为1,2,3时的一阶导，下方同理
    L=pabc/pall
  }else if(identical(xx,c(2,3,1))){
    L=pbca/pall
  }else if(identical(xx,c(1,3,2))){
    L=pacb/pall
  }else if(identical(xx,c(2,1,3))){
    L=pbac/pall
  }else if(identical(xx,c(3,1,2))){
    L=pcab/pall
  }else if(identical(xx,c(3,2,1))){
    L=pcba/pall
  }

  L
}

grad_theta_TIRT_3_rank <- function(theta,a,d,x=c(1,2,3)){

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
  pabc=p12*p13*p23
  pacb=p12*p13*p32
  pbac=p21*p23*p13
  pbca=p21*p23*p31
  pcab=p31*p32*p12
  pcba=p31*p32*p21
  pall=pabc+pacb+pbac+pbca+pcab+pcba
  if(identical(x,c(1,2,3))){       #顺序为1,2,3时的一阶导，下方同理
    g.i <- (ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p23+p12*p23*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))/pall-
      pabc/pall^2*((ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p23+p12*p23*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23*p31-p21*p23*ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2))+
                     (ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p32+p12*p32*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23*p13+p21*p23*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32*p12+p31*p32*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32*p21-p31*p32*ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)))
    g.j <- (-aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p23+p12*p13*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))/pall-
      pabc/pall^2*((-aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p23+p12*p13*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p23*p31+p21*p31*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (-aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p32-p12*p13*aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))+
                     (aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p23*p13+p21*p13*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (-aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)*p31*p12-p31*p32*aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2))+
                     (-aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)*p31*p21+p31*p32*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)))
    g.k <- (-p12*p23*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2)-p12*p13*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))/pall-
      pabc/pall^2*((-p12*p23*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2)-p12*p13*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (-p21*p31*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2)+p21*p23*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2))+
                     (-p12*p32*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2)+p12*p13*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))+
                     (-p21*p13*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2)-p21*p23*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (p32*p12*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)+p31*p12*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))+
                     (p32*p21*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)+p31*p21*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)))
    }else if(identical(x,c(2,3,1))){
    g.i <- (-ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23*p31-p21*p23*ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2))/pall-
      pbca/pall^2*((ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p23+p12*p23*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23*p31-p21*p23*ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2))+
                     (ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p32+p12*p32*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23*p13+p21*p23*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32*p12+p31*p32*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32*p21-p31*p32*ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)))
    g.j <- (aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p23*p31+p21*p31*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))/pall-
      pbca/pall^2*((-aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p23+p12*p13*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p23*p31+p21*p31*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (-aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p32-p12*p13*aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))+
                     (aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p23*p13+p21*p13*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (-aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)*p31*p12-p31*p32*aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2))+
                     (-aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)*p31*p21+p31*p32*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)))
    g.k <- (-p21*p31*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2)+p21*p23*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2))/pall-
      pbca/pall^2*((-p12*p23*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2)-p12*p13*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (-p21*p31*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2)+p21*p23*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2))+
                     (-p12*p32*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2)+p12*p13*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))+
                     (-p21*p13*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2)-p21*p23*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (p32*p12*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)+p31*p12*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))+
                     (p32*p21*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)+p31*p21*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)))
    }else if(identical(x,c(1,3,2))){
    g.i <- (ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p32+p12*p32*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))/pall-
      pacb/pall^2*((ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p23+p12*p23*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23*p31-p21*p23*ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2))+
                     (ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p32+p12*p32*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23*p13+p21*p23*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32*p12+p31*p32*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32*p21-p31*p32*ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)))
    g.j <- (-aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p32-p12*p13*aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))/pall-
      pacb/pall^2*((-aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p23+p12*p13*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p23*p31+p21*p31*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (-aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p32-p12*p13*aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))+
                     (aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p23*p13+p21*p13*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (-aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)*p31*p12-p31*p32*aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2))+
                     (-aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)*p31*p21+p31*p32*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)))
    g.k <- (-p12*p32*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2)+p12*p13*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))/pall-
      pacb/pall^2*((-p12*p23*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2)-p12*p13*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (-p21*p31*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2)+p21*p23*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2))+
                     (-p12*p32*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2)+p12*p13*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))+
                     (-p21*p13*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2)-p21*p23*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (p32*p12*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)+p31*p12*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))+
                     (p32*p21*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)+p31*p21*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)))
     }else if(identical(x,c(2,1,3))){
    g.i <- (-ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23*p13+p21*p23*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))/pall-
      pbac/pall^2*((ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p23+p12*p23*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23*p31-p21*p23*ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2))+
                     (ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p32+p12*p32*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23*p13+p21*p23*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32*p12+p31*p32*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32*p21-p31*p32*ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)))
    g.j <- (aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p23*p13+p21*p13*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))/pall-
      pbac/pall^2*((-aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p23+p12*p13*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p23*p31+p21*p31*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (-aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p32-p12*p13*aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))+
                     (aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p23*p13+p21*p13*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (-aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)*p31*p12-p31*p32*aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2))+
                     (-aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)*p31*p21+p31*p32*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)))
    g.k <- (-p21*p13*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2)-p21*p23*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))/pall-
      pbac/pall^2*((-p12*p23*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2)-p12*p13*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (-p21*p31*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2)+p21*p23*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2))+
                     (-p12*p32*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2)+p12*p13*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))+
                     (-p21*p13*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2)-p21*p23*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (p32*p12*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)+p31*p12*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))+
                     (p32*p21*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)+p31*p21*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)))
   }else if(identical(x,c(3,1,2))){
    g.i <- (-ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32*p12+p31*p32*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2))/pall-
      pcab/pall^2*((ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p23+p12*p23*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23*p31-p21*p23*ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2))+
                     (ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p32+p12*p32*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23*p13+p21*p23*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32*p12+p31*p32*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32*p21-p31*p32*ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)))
    g.j <- (-aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)*p31*p12-p31*p32*aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2))/pall-
      pcab/pall^2*((-aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p23+p12*p13*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p23*p31+p21*p31*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (-aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p32-p12*p13*aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))+
                     (aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p23*p13+p21*p13*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (-aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)*p31*p12-p31*p32*aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2))+
                     (-aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)*p31*p21+p31*p32*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)))
    g.k <- (p32*p12*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)+p31*p12*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))/pall-
      pcab/pall^2*((-p12*p23*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2)-p12*p13*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (-p21*p31*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2)+p21*p23*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2))+
                     (-p12*p32*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2)+p12*p13*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))+
                     (-p21*p13*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2)-p21*p23*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (p32*p12*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)+p31*p12*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))+
                     (p32*p21*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)+p31*p21*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)))
    }else if(identical(x,c(3,2,1))){
    g.i <- (-ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32*p21-p31*p32*ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2))/pall-
      pcba/pall^2*((ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p23+p12*p23*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23*p31-p21*p23*ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2))+
                     (ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p32+p12*p32*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23*p13+p21*p23*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32*p12+p31*p32*ai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2))+
                     (-ai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32*p21-p31*p32*ai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)))
    g.j <- (-aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)*p31*p21+p31*p32*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2))/pall-
      pcba/pall^2*((-aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p23+p12*p13*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p23*p31+p21*p31*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (-aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13*p32-p12*p13*aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))+
                     (aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p23*p13+p21*p13*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (-aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)*p31*p12-p31*p32*aj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2))+
                     (-aj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)*p31*p21+p31*p32*aj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)))
    g.k <- (p32*p21*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)+p31*p21*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))/pall-
      pcba/pall^2*((-p12*p23*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2)-p12*p13*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+
                     (-p21*p31*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2)+p21*p23*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2))+
                     (-p12*p32*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2)+p12*p13*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))+
                     (-p21*p13*ak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2)-p21*p23*ak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+
                     (p32*p12*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)+p31*p12*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))+
                     (p32*p21*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)+p31*p21*ak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)))
     }

  c(g.i,g.j,g.k)
}



Ib_theta_TIRT_3_rank <- function(theta,a,d){
  all.x <- list(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
  ret <- 0
  for(x in 1:length(all.x)){
    g <- grad_theta_TIRT_3_rank(theta,a,d,x=all.x[[x]])
    ret <- ret + tcrossprod(g,g)/px_2PL_3_rank(theta,a,d,xx=all.x[[x]])
  }
  ret
}


test.info_theta_TIRT_3_rank <- function(theta,item.par,BID){
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
    complete.info.b[dim.b,dim.b] <- info[[b]] <- Ib_theta_TIRT_3_rank(theta.b,a.b,d.b)
    I <- I+complete.info.b
    comp.info[[b]] <- complete.info.b
  }
  se <- sqrt(diag(solve(Reduce("+",comp.info))))
  list(comp.info=comp.info,se =se, info.b=info,I=I)
}
