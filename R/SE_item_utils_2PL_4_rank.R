px_item_2PL_4_rank <- function(theta,a,d,x=c(1,2,3,4)){ #a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)
  z <- exp(t(a*t(theta)-d))
  apply(z[,x[-4]], 1, prod)/(rowSums(z)*rowSums(z[,x[-1]])*rowSums(z[,x[-1:-2]]))
}
  grad_FDM_item_2PL_4_rank <- function(theta,a,d,x=c(1,2,3,4),h){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)

    g.ia <- (px_item_2PL_4_rank(theta,a+c(h,0,0,0),d,x)-px_item_2PL_4_rank(theta,a,d,x))/h
    g.ja <- (px_item_2PL_4_rank(theta,a+c(0,h,0,0),d,x)-px_item_2PL_4_rank(theta,a,d,x))/h
    g.ka <- (px_item_2PL_4_rank(theta,a+c(0,0,h,0),d,x)-px_item_2PL_4_rank(theta,a,d,x))/h
    g.la <- (px_item_2PL_4_rank(theta,a+c(0,0,0,h),d,x)-px_item_2PL_4_rank(theta,a,d,x))/h
    g.id <- (px_item_2PL_4_rank(theta,a,d+c(h,0,0,-h),x)-px_item_2PL_4_rank(theta,a,d,x))/h
    g.jd <- (px_item_2PL_4_rank(theta,a,d+c(0,h,0,-h),x)-px_item_2PL_4_rank(theta,a,d,x))/h
    g.kd <- (px_item_2PL_4_rank(theta,a,d+c(0,0,h,-h),x)-px_item_2PL_4_rank(theta,a,d,x))/h
    cbind(g.ia,g.ja,g.ka,g.la,g.id,g.jd,g.kd)
  }
  grad_CDM_item_2PL_4_rank <- function(theta,a,d,x=c(1,2,3,4),h){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)

    g.ia <- (px_item_2PL_4_rank(theta,a+c(h,0,0,0),d,x)-px_item_2PL_4_rank(theta,a-c(h,0,0,0),d,x))/2*h
    g.ja <- (px_item_2PL_4_rank(theta,a+c(0,h,0,0),d,x)-px_item_2PL_4_rank(theta,a-c(0,h,0,0),d,x))/2*h
    g.ka <- (px_item_2PL_4_rank(theta,a+c(0,0,h,0),d,x)-px_item_2PL_4_rank(theta,a-c(0,0,h,0),d,x))/2*h
    g.la <- (px_item_2PL_4_rank(theta,a+c(0,0,0,h),d,x)-px_item_2PL_4_rank(theta,a-c(0,0,0,h),d,x))/2*h
    g.id <- (px_item_2PL_4_rank(theta,a,d+c(h,0,0,-h),x)-px_item_2PL_4_rank(theta,a,d-c(h,0,0,-h),x))/2*h
    g.jd <- (px_item_2PL_4_rank(theta,a,d+c(0,h,0,-h),x)-px_item_2PL_4_rank(theta,a,d-c(0,h,0,-h),x))/2*h
    g.kd <- (px_item_2PL_4_rank(theta,a,d+c(0,0,h,-h),x)-px_item_2PL_4_rank(theta,a,d-c(0,0,h,-h),x))/2*h
    cbind(g.ia,g.ja,g.ka,g.la,g.id,g.jd,g.kd)
  }
  grad_RES_item_2PL_4_rank <- function(theta,a,d,x=c(1,2,3,4),h){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)

    g.ia <- (px_item_2PL_4_rank(theta,a-c(2*h,0,0,0),d,x)+8*px_item_2PL_4_rank(theta,a+c(h,0,0,0),d,x)
             -px_item_2PL_4_rank(theta,a+c(2*h,0,0,0),d,x)-8*px_item_2PL_4_rank(theta,a-c(h,0,0,0),d,x))/12*h
    g.ja <- (px_item_2PL_4_rank(theta,a-c(0,2*h,0,0),d,x)+8*px_item_2PL_4_rank(theta,a+c(0,h,0,0),d,x)
             -px_item_2PL_4_rank(theta,a+c(0,2*h,0,0),d,x)-8*px_item_2PL_4_rank(theta,a-c(0,h,0,0),d,x))/12*h
    g.ka <- (px_item_2PL_4_rank(theta,a-c(0,0,2*h,0),d,x)+8*px_item_2PL_4_rank(theta,a+c(0,0,h,0),d,x)
             -px_item_2PL_4_rank(theta,a+c(0,0,2*h,0),d,x)-8*px_item_2PL_4_rank(theta,a-c(0,0,h,0),d,x))/12*h
    g.la <- (px_item_2PL_4_rank(theta,a-c(0,0,0,2*h),d,x)+8*px_item_2PL_4_rank(theta,a+c(0,0,0,h),d,x)
             -px_item_2PL_4_rank(theta,a+c(0,0,0,2*h),d,x)-8*px_item_2PL_4_rank(theta,a-c(0,0,0,h),d,x))/12*h
    g.id <- (px_item_2PL_4_rank(theta,a,d-c(2*h,0,0,-2*h),x)+8*px_item_2PL_4_rank(theta,a,d+c(h,0,0,-h),x)
             -px_item_2PL_4_rank(theta,a,d+c(2*h,0,0,-2*h),x)-8*px_item_2PL_4_rank(theta,a,d-c(h,0,0,-h),x))/12*h
    g.jd <- (px_item_2PL_4_rank(theta,a,d-c(0,2*h,0,-2*h),x)+8*px_item_2PL_4_rank(theta,a,d+c(0,h,0,-h),x)
             -px_item_2PL_4_rank(theta,a,d+c(0,2*h,0,-2*h),x)-8*px_item_2PL_4_rank(theta,a,d-c(0,h,0,-h),x))/12*h
    g.kd <- (px_item_2PL_4_rank(theta,a,d-c(0,0,2*h,-2*h),x)+8*px_item_2PL_4_rank(theta,a,d+c(0,0,h,-h),x)
             -px_item_2PL_4_rank(theta,a,d+c(0,0,2*h,-2*h),x)-8*px_item_2PL_4_rank(theta,a,d-c(0,0,h,-h),x))/12*h
    cbind(g.ia,g.ja,g.ka,g.la,g.id,g.jd,g.kd)
  }

  grad_item_2PL_4_rank <- function(theta,a,d,x=c(1,2,3,4)){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)
    z <- exp(t(a*t(theta)-d))
    thetai <- theta[,1]
    thetaj <- theta[,2]
    thetak <- theta[,3]
    thetal <- theta[,4]
    # p.k.ijk <- p.k.ikj <- p.k.jki <- p.k.jik <- p.k.kij <- p.k.kji
    p.k <- z[,3]/rowSums(z)
    p.l <- z[,4]/rowSums(z)
    # p.j.ijk <- p.j.ikj <- p.j.jki <- p.j.jik <- p.j.kij <- p.j.kji
    p.j <- z[,2]/rowSums(z)
    # p.i.ijk <- p.i.ikj <- p.i.jki <- p.i.jik <- p.i.kij <- p.i.kji
    p.i <- z[,1]/rowSums(z)
    p.i.ijk <- p.i.ikj <- p.i.jik <- p.i.jki <- p.i.kij <- p.i.kji <- z[,1]/rowSums(z[,-4])
    p.j.ijk <- p.j.ikj <- p.j.jik <- p.j.jki <- p.j.kij <- p.j.kji <- z[,2]/rowSums(z[,-4])
    p.k.ijk <- p.k.ikj <- p.k.jik <- p.k.jki <- p.k.kij <- p.k.kji <- z[,3]/rowSums(z[,-4])
    p.i.ikl <- p.i.ilk <- p.i.kil <- p.i.kli <- p.i.lik <- p.i.lki <- z[,1]/rowSums(z[,-2])
    p.k.ikl <- p.k.ilk <- p.k.kil <- p.k.kli <- p.k.lik <- p.k.lki <- z[,3]/rowSums(z[,-2])
    p.l.ikl <- p.l.ilk <- p.l.kil <- p.l.kli <- p.l.lik <- p.l.lki <- z[,4]/rowSums(z[,-2])
    p.i.ijl <- p.i.ilj <- p.i.jil <- p.i.jli <- p.i.lij <- p.i.lji <- z[,1]/rowSums(z[,-3])
    p.j.ijl <- p.j.ilj <- p.j.jil <- p.j.jli <- p.j.lij <- p.j.lji <- z[,2]/rowSums(z[,-3])
    p.l.ijl <- p.l.ilj <- p.l.jil <- p.l.jli <- p.l.lij <- p.l.lji <- z[,4]/rowSums(z[,-3])
    p.j.jkl <- p.j.jlk <- p.j.klj <- p.j.kjl <- p.j.ljk <- p.j.lkj <- z[,2]/rowSums(z[,-1])
    p.k.jkl <- p.k.jlk <- p.k.klj <- p.k.kjl <- p.k.ljk <- p.k.lkj <- z[,3]/rowSums(z[,-1])
    p.l.jkl <- p.l.jlk <- p.l.klj <- p.l.kjl <- p.l.ljk <- p.l.lkj <- z[,4]/rowSums(z[,-1])

    p.k.kl <- p.k.lk <- z[,3]/rowSums(z[,-1:-2])
    p.l.kl <- p.l.lk <- z[,4]/rowSums(z[,-1:-2])
    p.i.ij <- p.i.ji <- z[,1]/rowSums(z[,-3:-4])
    p.j.ij <- p.j.ji <- z[,2]/rowSums(z[,-3:-4])
    p.i.ik <- p.i.ki <- z[,1]/rowSums(z[,c(-2,-4)])
    p.k.ik <- p.k.ki <- z[,3]/rowSums(z[,c(-2,-4)])
    p.i.il <- p.i.li <- z[,1]/rowSums(z[,-2:-3])
    p.l.il <- p.l.li <- z[,4]/rowSums(z[,-2:-3])
    p.j.jk <- p.j.kj <- z[,2]/rowSums(z[,c(-1,-4)])
    p.k.jk <- p.k.kj <- z[,3]/rowSums(z[,c(-1,-4)])
    p.j.jl <- p.j.lj <- z[,2]/rowSums(z[,c(-1,-3)])
    p.l.jl <- p.l.lj <- z[,4]/rowSums(z[,c(-1,-3)])
    if(identical(x,c(1,2,3,4))){
      g.ia <- thetai*p.i*(1-p.i)*p.j.jkl*p.k.kl
      g.ja <- thetaj*p.i*p.j.jkl*p.k.kl*(1-p.j-p.j.jkl)
      g.ka <- thetak*p.i*p.j.jkl*p.k.kl*(p.l.kl-p.k.jkl-p.k)
      g.la <- -thetal*p.i*p.j.jkl*p.k.kl*(p.l+p.l.jkl+p.l.kl)
      g.id <- -p.i*p.j.jkl*p.k.kl*(p.j+p.k+2*p.l+p.l.jkl+p.l.kl)
      g.jd <- p.i*p.j.jkl*p.k.kl*(p.j-p.l-p.k.jkl-2*p.l.jkl-p.l.kl)
      g.kd <- p.i*p.j.jkl*p.k.kl*(p.k-p.l+p.k.jkl-p.l.jkl-2*p.l.kl)
    }else if(identical(x,c(1,2,4,3))){
      g.ia <- thetai*p.i*(1-p.i)*p.j.jlk*p.l.lk
      g.ja <- thetaj*p.i*p.j.jlk*p.l.lk*(1-p.j-p.j.jlk)
      g.la <- thetal*p.i*p.j.jlk*p.l.lk*(p.k.lk-p.l.jlk-p.l)
      g.ka <- -thetak*p.i*p.j.jlk*p.l.lk*(p.k+p.k.jlk+p.k.lk)
      g.id <- -p.i*p.j.jlk*p.l.lk*(p.j+p.k+2*p.l+p.l.jlk-p.k.lk)
      g.jd <- p.i*p.j.jlk*p.l.lk*(p.j-p.l-2*p.l.jlk-p.k.jkl+p.k.lk)
      g.kd <- p.i*p.j.jlk*p.l.lk*(p.k-p.l+p.k.jlk-p.l.jlk+2*p.k.lk)
    }else if(identical(x,c(1,3,2,4))){
      g.ia <- thetai*p.i*(1-p.i)*p.k.kjl*p.j.jl
      g.ka <- thetak*p.i*p.k.kjl*p.j.jl*(1-p.k-p.k.kjl)
      g.ja <- thetaj*p.i*p.k.kjl*p.j.jl*(p.l.jl-p.j.kjl-p.j)
      g.la <- -thetal*p.i*p.k.kjl*p.j.jl*(p.l+p.l.kjl+p.l.jl)
      g.id <- -p.i*p.j.jkl*p.k.kl*(p.j+p.k+2*p.l+p.l.jkl+p.l.kl)
      g.jd <- p.i*p.k.jkl*p.j.jl*(p.k-p.l-p.j.jkl-2*p.l.jkl-p.l.jl)
      g.kd <- p.i*p.k.jkl*p.j.jl*(p.j-p.l+p.j.jkl-p.l.jkl-2*p.l.jl)
    }else if(identical(x,c(1,3,4,2))){
      g.ia <- thetai*p.i*(1-p.i)*p.k.klj*p.l.lj
      g.ka <- thetak*p.i*p.k.klj*p.l.lj*(1-p.k-p.k.klj)
      g.la <- thetal*p.i*p.k.klj*p.l.lj*(p.j.lj-p.l.klj-p.l)
      g.ja <- -thetaj*p.i*p.k.klj*p.l.lj*(p.j+p.j.klj+p.j.lj)
      g.id <- -p.i*p.k.klj*p.l.lj*(p.k+p.j+2*p.l+p.l.klj-p.j.lj)
      g.kd <- p.i*p.k.klj*p.l.lj*(p.k-p.l-2*p.l.klj-p.j.kjl+p.j.lj)
      g.jd <- p.i*p.k.klj*p.l.lj*(p.j-p.l+p.j.klj-p.l.klj+2*p.j.lj)
    }else if(identical(x,c(1,4,2,3))){
      g.ia <- thetai*p.i*(1-p.i)*p.l.ljk*p.j.jk
      g.la <- thetal*p.i*p.l.ljk*p.j.jk*(1-p.l-p.l.ljk)
      g.ja <- thetaj*p.i*p.l.ljk*p.j.jk*(p.k.jk-p.j.ljk-p.j)
      g.ka <- -thetak*p.i*p.l.ljk*p.j.jk*(p.k+p.k.ljk+p.k.jk)
      g.id <- -p.i*p.l.ljk*p.j.jk*(2*p.l+p.j+p.k-p.j.ljk-p.k.ljk)
      g.jd <- p.i*p.l.ljk*p.j.jk*(p.j-p.l-p.k.jk+2*p.j.ljk+p.k.ljk)
      g.kd <- p.i*p.l.ljk*p.j.jk*(p.k-p.l+2*p.k.ljk+p.j.ljk+p.k.jk)
    }else if(identical(x,c(1,4,3,2))){
      g.ia <- thetai*p.i*(1-p.i)*p.l.lkj*p.k.kj
      g.la <- thetal*p.i*p.l.lkj*p.k.kj*(1-p.l-p.l.lkj)
      g.ka <- thetak*p.i*p.l.lkj*p.k.kj*(p.j.kj-p.k.lkj-p.k)
      g.ja <- -thetaj*p.i*p.l.lkj*p.k.kj*(p.j+p.j.lkj+p.j.kj)
      g.id <- -p.i*p.l.lkj*p.k.kj*(2*p.l+p.k+p.j-p.k.lkj-p.j.lkj)
      g.kd <- p.i*p.l.lkj*p.k.kj*(p.k-p.l-p.j.kj+2*p.k.lkj+p.j.lkj)
      g.jd <- p.i*p.l.lkj*p.k.kj*(p.j-p.l+2*p.j.lkj+p.k.lkj+p.j.kj)
    }else if(identical(x,c(2,1,3,4))){
      g.ja <- thetaj*p.j*(1-p.j)*p.i.ikl*p.k.kl
      g.ia <- thetai*p.j*p.i.ikl*p.k.kl*(1-p.i-p.i.ikl)
      g.ka <- thetak*p.j*p.i.ikl*p.k.kl*(p.l.kl-p.k.ikl-p.k)
      g.la <- -thetal*p.j*p.i.ikl*p.k.kl*(p.l+p.l.ikl+p.l.kl)
      g.id <- -p.j*p.i.ikl*p.k.kl*(p.i+p.k+2*p.l+p.l.ikl+p.l.kl)
      g.jd <- p.j*p.i.ikl*p.k.kl*(p.i-p.l-p.k.ikl-2*p.l.ikl-p.l.kl)
      g.kd <- p.i*p.j.jkl*p.k.kl*(p.k-p.l+p.k.jkl-p.l.jkl-2*p.l.kl)
    }else if(identical(x,c(2,1,4,3))){
      g.ja <- thetaj*p.j*(1-p.j)*p.i.ilk*p.l.lk
      g.ia <- thetai*p.j*p.i.ilk*p.l.lk*(1-p.i-p.i.ilk)
      g.la <- thetal*p.j*p.i.ilk*p.l.lk*(p.k.lk-p.l.ilk-p.l)
      g.ka <- -thetak*p.j*p.i.ilk*p.l.lk*(p.k+p.k.ilk+p.k.lk)
      g.jd <- -p.j*p.i.ikl*p.k.kl*(p.i+p.k+2*p.l+p.l.ikl+p.l.kl)
      g.id <- p.j*p.i.ikl*p.k.kl*(p.i-p.l-p.k.ikl-2*p.l.ikl-p.l.kl)
      g.kd <- p.j*p.i.ikl*p.k.kl*(p.k-p.l+p.k.ikl-p.l.ikl-2*p.l.kl)
    }else if(identical(x,c(2,3,1,4))){
      g.ja <- thetaj*p.j*(1-p.j)*p.k.kil*p.i.il
      g.ka <- thetak*p.j*p.k.kil*p.i.il*(1-p.k-p.k.kil)
      g.ia <- thetai*p.j*p.k.kil*p.i.il*(p.l.il-p.i.kil-p.i)
      g.la <- -thetal*p.j*p.k.kil*p.i.il*(p.l+p.l.kil+p.l.il)
      g.jd <- -p.j*p.k.kil*p.i.il*(p.k+p.i+2*p.l+p.l.kil+p.l.il)
      g.kd <- p.j*p.k.kil*p.i.il*(p.k-p.l-p.i.kil-2*p.l.kil-p.l.il)
      g.id <- p.j*p.k.kil*p.i.il*(p.i-p.l+p.i.kil-p.l.kil-2*p.l.il)
    }else if(identical(x,c(2,3,4,1))){
      g.ja <- thetaj*p.j*(1-p.j)*p.k.kli*p.l.li
      g.ka <- thetak*p.j*p.k.kli*p.l.li*(1-p.k-p.k.kli)
      g.la <- thetal*p.j*p.k.kli*p.l.li*(p.i.li-p.l.kli-p.l)
      g.ia <- -thetai*p.j*p.k.kli*p.l.li*(p.i+p.i.kli+p.i.li)
      g.jd <- -p.j*p.k.kil*p.i.il*(p.k+p.i+2*p.l+p.l.kil+p.l.il)
      g.kd <- p.j*p.k.kil*p.i.il*(p.k-p.l-p.i.kil-2*p.l.kil-p.l.il)
      g.id <- p.j*p.k.kil*p.i.il*(p.i-p.l+p.i.kil-p.l.kil-2*p.l.il)
    }else if(identical(x,c(2,4,1,3))){
      g.ja <- thetaj*p.j*(1-p.j)*p.l.lik*p.i.ik
      g.la <- thetal*p.j*p.l.lik*p.i.ik*(1-p.l-p.l.lik)
      g.ia <- thetai*p.j*p.l.lik*p.i.ik*(p.k.ik-p.i.lik-p.i)
      g.ka <- -thetak*p.j*p.l.lik*p.i.ik*(p.k+p.k.lik+p.k.ik)
      g.jd <- -p.j*p.l.lik*p.i.ik*(2*p.l+p.i+p.k-p.i.lik-p.k.lik)
      g.id <- p.j*p.l.lik*p.i.ik*(p.i-p.l-p.k.ik+2*p.i.lik+p.k.lik)
      g.kd <- p.j*p.l.lik*p.i.ik*(p.k-p.l+2*p.k.lik+p.i.lik+p.k.ik)
    }else if(identical(x,c(2,4,3,1))){
      g.ja <- thetaj*p.j*(1-p.j)*p.l.lki*p.k.ki
      g.la <- thetal*p.j*p.l.lki*p.k.ki*(1-p.l-p.l.lki)
      g.ka <- thetak*p.j*p.l.lki*p.k.ki*(p.i.ki-p.k.lki-p.k)
      g.ia <- -thetai*p.j*p.l.lki*p.k.ki*(p.i+p.i.lki+p.i.ki)
      g.jd <- -p.j*p.l.lki*p.k.ki*(2*p.l+p.k+p.i-p.k.lki-p.i.lki)
      g.kd <- p.j*p.l.lki*p.k.ki*(p.k-p.l-p.i.ki+2*p.k.lki+p.i.lki)
      g.id <- p.j*p.l.lki*p.k.ki*(p.i-p.l+2*p.i.lki+p.k.lki+p.i.ki)
    }else if(identical(x,c(3,1,2,4))){
      g.ka <- thetak*p.k*(1-p.k)*p.i.ijl*p.j.jl
      g.ia <- thetai*p.k*p.i.ijl*p.j.jl*(1-p.i-p.i.ijl)
      g.ja <- thetaj*p.k*p.i.ijl*p.j.jl*(p.l.jl-p.j.ijl-p.j)
      g.la <- -thetal*p.k*p.i.ijl*p.j.jl*(p.l+p.l.ijl+p.l.jl)
      g.kd <- -p.k*p.i.ijl*p.j.jl*(p.i+p.j+2*p.l+p.l.ijl+p.l.jl)
      g.id <- p.k*p.i.ijl*p.j.jl*(p.i-p.l-p.j.ijl-2*p.l.ijl-p.l.jl)
      g.jd <- p.k*p.i.ijl*p.j.jl*(p.j-p.l+p.j.ijl-p.l.ijl-2*p.l.jl)
    }else if(identical(x,c(3,1,4,2))){
      g.ka <- thetak*p.k*(1-p.k)*p.i.ilj*p.l.lj
      g.ia <- thetai*p.k*p.i.ilj*p.l.lj*(1-p.i-p.i.ilj)
      g.la <- thetal*p.k*p.i.ilj*p.l.lj*(p.j.lj-p.l.ilj-p.l)
      g.ja <- -thetaj*p.k*p.i.ilj*p.l.lj*(p.j+p.j.ilj+p.j.lj)
      g.kd <- -p.k*p.i.ijl*p.j.jl*(p.i+p.j+2*p.l+p.l.ijl+p.l.jl)
      g.id <- p.k*p.i.ijl*p.j.jl*(p.i-p.l-p.j.ijl-2*p.l.ijl-p.l.jl)
      g.jd <- p.k*p.i.ijl*p.j.jl*(p.j-p.l+p.j.ijl-p.l.ijl-2*p.l.jl)
    }else if(identical(x,c(3,2,1,4))){
      g.ka <- thetak*p.k*(1-p.k)*p.j.jil*p.i.il
      g.ja <- thetaj*p.k*p.j.jil*p.i.il*(1-p.j-p.j.jil)
      g.ia <- thetai*p.k*p.j.jil*p.i.il*(p.l.il-p.i.jil-p.i)
      g.la <- -thetal*p.k*p.j.jil*p.i.il*(p.l+p.l.jil+p.l.il)
      g.kd <- -p.k*p.j.jil*p.i.il*(p.j+p.i+2*p.l+p.l.jil+p.l.il)
      g.jd <- p.k*p.j.jil*p.i.il*(p.j-p.l-p.i.jil-2*p.l.jil-p.l.il)
      g.id <- p.k*p.j.jil*p.i.il*(p.i-p.l+p.i.jil-p.l.jil-2*p.l.il)
    }else if(identical(x,c(3,2,4,1))){
      g.ka <- thetak*p.k*(1-p.k)*p.j.jli*p.l.li
      g.ja <- thetaj*p.k*p.j.jli*p.l.li*(1-p.j-p.j.jli)
      g.la <- thetal*p.k*p.j.jli*p.l.li*(p.i.li-p.l.jli-p.l)
      g.ia <- -thetai*p.k*p.j.jli*p.l.li*(p.i+p.i.jli+p.i.li)
      g.kd <- -p.k*p.j.jil*p.i.il*(p.j+p.i+2*p.l+p.l.jil+p.l.il)
      g.jd <- p.k*p.j.jil*p.i.il*(p.j-p.l-p.i.jil-2*p.l.jil-p.l.il)
      g.id <- p.k*p.j.jil*p.i.il*(p.i-p.l+p.i.jil-p.l.jil-2*p.l.il)
    }else if(identical(x,c(3,4,1,2))){
      g.ka <- thetak*p.k*(1-p.k)*p.l.lij*p.i.ij
      g.la <- thetal*p.k*p.l.lij*p.i.ij*(1-p.l-p.l.lij)
      g.ia <- thetai*p.k*p.l.lij*p.i.ij*(p.j.ij-p.i.lij-p.i)
      g.ja <- -thetaj*p.k*p.l.lij*p.i.ij*(p.j+p.j.lij+p.j.ij)
      g.kd <- -p.k*p.l.lij*p.i.ij*(2*p.l+p.i+p.j-p.i.lij-p.j.lij)
      g.id <- p.k*p.l.lij*p.i.ij*(p.i-p.l-p.j.ij+2*p.i.lij+p.j.lij)
      g.jd <- p.k*p.l.lij*p.i.ij*(p.j-p.l+2*p.j.lij+p.i.lij+p.j.ij)
    }else if(identical(x,c(3,4,2,1))){
      g.ka <- thetak*p.k*(1-p.k)*p.l.lji*p.j.ji
      g.la <- thetal*p.k*p.l.lji*p.j.ji*(1-p.l-p.l.lji)
      g.ja <- thetaj*p.k*p.l.lji*p.j.ji*(p.i.ji-p.j.lji-p.j)
      g.ia <- -thetai*p.k*p.l.lji*p.j.ji*(p.i+p.i.lji+p.i.ji)
      g.kd <- -p.k*p.l.lji*p.j.ji*(2*p.l+p.j+p.i-p.j.lji-p.i.lji)
      g.jd <- p.k*p.l.lji*p.j.ji*(p.j-p.l-p.i.ji+2*p.j.lji+p.i.lji)
      g.id <- p.k*p.l.lji*p.j.ji*(p.i-p.l+2*p.i.lji+p.j.lji+p.i.ji)
    }else if(identical(x,c(4,1,2,3))){
      g.la <- thetal*p.l*(1-p.l)*p.i.ijk*p.j.jk
      g.ia <- thetai*p.l*p.i.ijk*p.j.jk*(1-p.i-p.i.ijk)
      g.ja <- thetaj*p.l*p.i.ijk*p.j.jk*(p.k.jk-p.j.ijk-p.j)
      g.ka <- -thetak*p.l*p.i.ijk*p.j.jk*(p.k+p.k.ijk+p.k.jk)
      g.id <- p.l*p.i.ijk*p.j.jk*(2*p.i+p.j+p.k-p.j.ijk-p.k.ijk)
      g.jd <- p.l*p.i.ijk*p.j.jk*(p.i+2*p.j+p.k+p.j.ijk-p.k.jk)
      g.kd <- p.l*p.i.ijk*p.j.jk*(p.i+p.j+2*p.k+p.k.ijk+p.k.jk)
    }else if(identical(x,c(4,1,3,2))){
      g.la <- thetal*p.l*(1-p.l)*p.i.ikj*p.k.kj
      g.ia <- thetai*p.l*p.i.ikj*p.k.kj*(1-p.i-p.i.ikj)
      g.ka <- thetak*p.l*p.i.ikj*p.k.kj*(p.j.kj-p.k.ikj-p.k)
      g.ja <- -thetaj*p.l*p.i.ikj*p.k.kj*(p.j+p.j.ikj+p.j.kj)
      g.id <- p.l*p.i.ikj*p.k.kj*(2*p.i+p.k+p.j-p.k.ikj-p.j.ikj)
      g.kd <- p.l*p.i.ikj*p.k.kj*(p.i+2*p.k+p.j+p.k.ikj-p.j.kj)
      g.jd <- p.l*p.i.ikj*p.k.kj*(p.i+p.k+2*p.j+p.j.ikj+p.j.kj)
    }else if(identical(x,c(4,2,1,3))){
      g.la <- thetal*p.l*(1-p.l)*p.j.jik*p.i.ik
      g.ja <- thetaj*p.l*p.j.jik*p.i.ik*(1-p.j-p.j.jik)
      g.ia <- thetai*p.l*p.j.jik*p.i.ik*(p.k.ik-p.i.jik-p.i)
      g.ka <- -thetak*p.l*p.j.jik*p.i.ik*(p.k+p.k.jik+p.k.ik)
      g.jd <- p.l*p.j.jik*p.i.ik*(2*p.j+p.i+p.k-p.i.jik-p.k.jik)
      g.id <- p.l*p.j.jik*p.i.ik*(p.j+2*p.i+p.k+p.i.jik-p.k.ik)
      g.kd <- p.l*p.j.jik*p.i.ik*(p.j+p.i+2*p.k+p.k.jik+p.k.ik)
    }else if(identical(x,c(4,2,3,1))){
      g.la <- thetal*p.l*(1-p.l)*p.j.jki*p.k.ki
      g.ja <- thetaj*p.l*p.j.jki*p.k.ki*(1-p.j-p.j.jki)
      g.ka <- thetak*p.l*p.j.jki*p.k.ki*(p.i.ki-p.k.jki-p.k)
      g.ia <- -thetai*p.l*p.j.jki*p.k.ki*(p.i+p.i.jki+p.i.ki)
      g.jd <- p.l*p.j.jki*p.k.ki*(2*p.j+p.k+p.i-p.k.jki-p.i.jki)
      g.kd <- p.l*p.j.jki*p.k.ki*(p.j+2*p.k+p.i+p.k.jki-p.i.ki)
      g.id <- p.l*p.j.jki*p.k.ki*(p.j+p.k+2*p.i+p.i.jki+p.i.ki)
    }else if(identical(x,c(4,3,1,2))){
      g.la <- thetal*p.l*(1-p.l)*p.k.kij*p.i.ij
      g.ka <- thetak*p.l*p.k.kij*p.i.ij*(1-p.k-p.k.kij)
      g.ia <- thetai*p.l*p.k.kij*p.i.ij*(p.j.ij-p.i.kij-p.i)
      g.ja <- -thetaj*p.l*p.k.kij*p.i.ij*(p.j+p.j.kij+p.j.ij)
      g.kd <- p.l*p.k.kij*p.i.ij*(2*p.k+p.i+p.j-p.i.kij-p.j.kij)
      g.id <- p.l*p.k.kij*p.i.ij*(p.k+2*p.i+p.j+p.i.kij-p.j.ij)
      g.jd <- p.l*p.k.kij*p.i.ij*(p.k+p.i+2*p.j+p.j.kij+p.j.ij)
    }else if(identical(x,c(4,3,2,1))){
      g.la <- thetal*p.l*(1-p.l)*p.k.kji*p.j.ji
      g.ka <- thetak*p.l*p.k.kji*p.j.ji*(1-p.k-p.k.kji)
      g.ja <- thetaj*p.l*p.k.kji*p.j.ji*(p.i.ji-p.j.kji-p.j)
      g.ia <- -thetai*p.l*p.k.kji*p.j.ji*(p.i+p.i.kji+p.i.ji)
      g.kd <- p.l*p.k.kji*p.j.ji*(2*p.k+p.j+p.i-p.j.kji-p.i.kji)
      g.jd <- p.l*p.k.kji*p.j.ji*(p.k+2*p.j+p.i+p.j.kji-p.i.ji)
      g.id <- p.l*p.k.kji*p.j.ji*(p.k+p.j+2*p.i+p.i.kji+p.i.ji)
    }

    cbind(g.ia,g.ja,g.ka,g.la,g.id,g.jd,g.kd)
  }



  FDM_item_2PL_4_rank <- function(theta,a,d,Y,h=1e-5){
    all.x <- list(c(1,2,3,4),c(1,2,4,3),c(1,3,2,4),c(1,3,4,2),c(1,4,2,3),c(1,4,3,2),
                  c(2,1,3,4),c(2,1,4,3),c(2,3,1,4),c(2,3,4,1),c(2,4,1,3),c(2,4,3,1),
                  c(3,1,2,4),c(3,1,4,2),c(3,2,1,4),c(3,2,4,1),c(3,4,1,2),c(3,4,2,1),
                  c(4,1,2,3),c(4,1,3,2),c(4,2,1,3),c(4,2,3,1),c(4,3,1,2),c(4,3,2,1))
    g <- grad_FDM_item_2PL_4_rank(theta,a,d,x=all.x[[Y]],h=h)/px_item_2PL_4_rank(theta,a,d,x=all.x[[Y]])
    g
  }
  RES_item_2PL_4_rank <- function(theta,a,d,Y,h=1){
    all.x <- list(c(1,2,3,4),c(1,2,4,3),c(1,3,2,4),c(1,3,4,2),c(1,4,2,3),c(1,4,3,2),
                  c(2,1,3,4),c(2,1,4,3),c(2,3,1,4),c(2,3,4,1),c(2,4,1,3),c(2,4,3,1),
                  c(3,1,2,4),c(3,1,4,2),c(3,2,1,4),c(3,2,4,1),c(3,4,1,2),c(3,4,2,1),
                  c(4,1,2,3),c(4,1,3,2),c(4,2,1,3),c(4,2,3,1),c(4,3,1,2),c(4,3,2,1))
    g <- grad_RES_item_2PL_4_rank(theta,a,d,x=all.x[[Y]],h=h)/px_item_2PL_4_rank(theta,a,d,x=all.x[[Y]])
    g
  }
  CDM_item_2PL_4_rank <- function(theta,a,d,Y,h=1){
    all.x <- list(c(1,2,3,4),c(1,2,4,3),c(1,3,2,4),c(1,3,4,2),c(1,4,2,3),c(1,4,3,2),
                  c(2,1,3,4),c(2,1,4,3),c(2,3,1,4),c(2,3,4,1),c(2,4,1,3),c(2,4,3,1),
                  c(3,1,2,4),c(3,1,4,2),c(3,2,1,4),c(3,2,4,1),c(3,4,1,2),c(3,4,2,1),
                  c(4,1,2,3),c(4,1,3,2),c(4,2,1,3),c(4,2,3,1),c(4,3,1,2),c(4,3,2,1))
    g <- grad_CDM_item_2PL_4_rank(theta,a,d,x=all.x[[Y]],h=h)/px_item_2PL_4_rank(theta,a,d,x=all.x[[Y]])
    g
  }
  XPD_item_2PL_4_rank <- function(theta,a,d,Y){
    all.x <- list(c(1,2,3,4),c(1,2,4,3),c(1,3,2,4),c(1,3,4,2),c(1,4,2,3),c(1,4,3,2),
                  c(2,1,3,4),c(2,1,4,3),c(2,3,1,4),c(2,3,4,1),c(2,4,1,3),c(2,4,3,1),
                  c(3,1,2,4),c(3,1,4,2),c(3,2,1,4),c(3,2,4,1),c(3,4,1,2),c(3,4,2,1),
                  c(4,1,2,3),c(4,1,3,2),c(4,2,1,3),c(4,2,3,1),c(4,3,1,2),c(4,3,2,1))
    g <- grad_item_2PL_4_rank(theta,a,d,x=all.x[[Y]])/px_item_2PL_4_rank(theta,a,d,x=all.x[[Y]])
    g
  }


  itemI_2PL_4_rank <- function(theta,item.par,BID,Y,I,h=NULL){
    D <- max(BID$Block)*(max(BID$Item)*2-1) #number of dims
    itemI <- matrix(0,nrow=nrow(theta),ncol=D)
    for(b in unique(BID$Block)){
      items <- which(BID$Block==b)
      item <- c((b-1)*7+1,(b-1)*7+2,(b-1)*7+3,(b-1)*7+4,(b-1)*7+5,(b-1)*7+6,(b-1)*7+7)
      dim.b <- BID$Dim[items]
      a.b <- item.par[items,1]
      d.b <- item.par[items,2]
      theta.b <- theta[,dim.b]
      Y.b <- Y[b]
      if(I=="XPD"){
        itemI[,item] <- XPD_item_2PL_4_rank(theta.b,a.b,d.b,Y.b)
      }else if(I=="FDM"){
        if(is.null(h)){
          itemI[,item] <- FDM_item_2PL_4_rank(theta.b,a.b,d.b,Y.b)
        }else{
          itemI[,item] <- FDM_item_2PL_4_rank(theta.b,a.b,d.b,Y.b,h=h)
        }
      }else if(I=="CDM"){
        if(is.null(h)){
          itemI[,item] <- CDM_item_2PL_4_rank(theta.b,a.b,d.b,Y.b)
        }else{
          itemI[,item] <- CDM_item_2PL_4_rank(theta.b,a.b,d.b,Y.b,h=h)
        }
      }else if(I=="RES"){
        if(is.null(h)){
          itemI[,item] <- RES_item_2PL_4_rank(theta.b,a.b,d.b,Y.b)
        }else{
          itemI[,item] <- RES_item_2PL_4_rank(theta.b,a.b,d.b,Y.b,h=h)
        }
      }
    }
    itemI
  }


  itemSE_Fisher_2PL_4_rank <- function(x,BID,Y,IN,q=10,h=NULL){
    s1=Sys.time()
    nCPUcores = detectCores()
    cl = makeCluster(nCPUcores-1)
    registerDoParallel(cl)
    D <- max(BID$Dim)
    N <- nrow(Y)  # number of examinees
    J <- ncol(Y)  # number of items
    X <- seq(-6, 6, length.out = q)  # quadrature points
    AX <- dnorm(X, 0, 1) / sum(dnorm(X, 0, 1))  # quadrature approximation
    X_matrix=as.matrix(expand.grid(rep(list(X), D)))
    AX_matrix=as.matrix(expand.grid(rep(list(AX), D)))
    AX_p=apply(AX_matrix, 1, prod)

    #估计值或真值
    a=x$a
    d=x$d

    resorder=as.matrix(expand.grid(rep(list(1:K), J)))

    paiH1 <- matrix(0, nrow = nrow(resorder), ncol = nrow(resorder))
    p1_5=matrix(0,nrow=q^D,ncol=J)

    sinI_vector=matrix(0, nrow = IN*J, ncol = nrow(resorder))
    all.x <- list(c(1,2,3,4),c(1,2,4,3),c(1,3,2,4),c(1,3,4,2),c(1,4,2,3),c(1,4,3,2),
                  c(2,1,3,4),c(2,1,4,3),c(2,3,1,4),c(2,3,4,1),c(2,4,1,3),c(2,4,3,1),
                  c(3,1,2,4),c(3,1,4,2),c(3,2,1,4),c(3,2,4,1),c(3,4,1,2),c(3,4,2,1),
                  c(4,1,2,3),c(4,1,3,2),c(4,2,1,3),c(4,2,3,1),c(4,3,1,2),c(4,3,2,1))
    for (h in 1:nrow(resorder)) {
      newresorder1=matrix(rep(as.matrix(resorder)[h,], each = q^D), nrow = q^D, ncol = J, byrow = F)
      for(i in 1:q^D){
        for(j in 1:J){
          p1_5[i,j]=px_item_2PL_4_rank(X_matrix[i,],a[,j],d[,j],x=all.x[[newresorder1[i,j]]])
        }
      }
      paiH1[h,h]=sum(apply(p1_5, 1, prod)*AX_p)
      post=apply(p1_5, 1, prod)*AX_p/paiH1[h,h]

      sinI_vector[,h]=foreach(i = seq_len(nrow(X_matrix)),.inorder=FALSE,.combine='+',.export = c("itemI_2PL_4_rank",                                                                                                                                                                                    "grad_item_2PL_4_rank",
                                                                                                  "px_item_2PL_4_rank",
                                                                                                  "XPD_item_2PL_4_rank"))%dopar%{
                                                                                                    itemI_2PL_4_rank(X_matrix[i,],data.frame(a=c(x$a),d=c(x$d)),BID,resorder[h,],I='XPD',h=NULL)*post[i]
                                                                                                  }
    }
    I=sinI_vector%*%paiH1%*%t(sinI_vector)*N
    se.a=NULL
    se.d=NULL
    se.sigm=NULL
    tryCatch({
      se <- matrix(sqrt(diag(solve(I))),nrow=max(BID$Block),ncol=IN,byrow = TRUE)
      se.a <- se[,1:4]
      se.d <- se[,5:7]
      se.d <- cbind(se.d,sqrt(se.d[,1]^2+se.d[,2]^2+se.d[,3]^2))
    }, error = function(e) {
      cat("The information matrix is irreversible!")
    })
    s2=Sys.time()
    list(SE.a=se.a,SE.d=se.d,SE.sigm=se.sigm,I=I,time=s2-s1)
  }



  itemSE_Simple_2PL_4_rank <- function(x,BID,Y,SE='XPD',IN,hess=NULL,h1=NULL){
    s1=Sys.time()
    D <- max(BID$Dim)
    N <- nrow(Y)  # number of examinees
    J <- ncol(Y)  # number of items
    #估计值或真值
    a=x$a
    d=x$d
    thetaex=do.call(rbind, x$thetalist)


    paiH1 <- matrix(0, nrow = nrow(Y), ncol = nrow(Y))
    diag(paiH1)=1
    sinI_vector=matrix(0, nrow = IN*J, ncol = nrow(Y))
    sinI_vector1=matrix(0, nrow = IN*J, ncol = nrow(Y))
    IgradPP0=matrix(0, nrow = IN*J, ncol = IN*J)
    for (h in 1:nrow(Y)) {{
      if(h %% 200==0){
        cat(h,'/',nrow(Y),'\n')
      }
      step_value =nrow(x$thetalist[[1]])
      step_length=length(x$thetalist)
      # 生成新向量
      new_vector <- sapply(which(apply(Y, 1, function(x) all(x == as.matrix(Y[h,])))), function(x) seq(from = x, by = step_value, length.out = step_length))
      thetaexchoose=thetaex[c(new_vector),]


      if(SE=='Louis'|SE=='Sandwich'){
        g=itemI_2PL_4_rank(thetaexchoose,data.frame(a=c(x$a),d=c(x$d)),BID,as.matrix(Y[h,]),I='XPD',h=NULL)
        # 初始化结果矩阵
        result_sum <- matrix(0, nrow = ncol(g), ncol = ncol(g))
        # 对每一行计算外积并累加
        for (i in 1:nrow(g)) {
          vec <- g[i, ] # 提取第 i 行
          result_sum <- result_sum + outer(vec, vec) # 计算外积并累加
        }
        sinI_vector1 <- result_sum / nrow(g)

        sinI_vector[,h]=colMeans(g)
        IgradPP0=IgradPP0+sinI_vector1
      }else{
        sinI_vector[,h]=colMeans(itemI_2PL_4_rank(thetaexchoose,data.frame(a=c(x$a),d=c(x$d)),BID,as.matrix(Y[h,]),I=SE,h=h1))

      }
    }}

    IXPD=sinI_vector%*%paiH1%*%t(sinI_vector)
    I=IXPD
    if(SE=='Louis'|SE=='Sandwich'){
      ILouis=IXPD+hess-IgradPP0
      I=ILouis
      if(SE=='Sandwich'){
        Isandwich=ILouis%*%solve(IXPD)%*%ILouis
        I=Isandwich
      }

    }
    se.a=NULL
    se.d=NULL
    se.sigm=NULL
    tryCatch({
      se <- matrix(sqrt((diag(solve(I)))),nrow=max(BID$Block),ncol=IN,byrow = TRUE)
      se.a <- se[,1:4]
      se.d <- se[,5:7]
      se.d <- cbind(se.d,sqrt(se.d[,1]^2+se.d[,2]^2+se.d[,3]^2))
    }, error = function(e) {
      cat("The information matrix is irreversible!")
    })
    s2=Sys.time()
    list(SE.a=se.a,SE.d=se.d,SE.sigm=se.sigm,I=I,time=s2-s1)
  }





