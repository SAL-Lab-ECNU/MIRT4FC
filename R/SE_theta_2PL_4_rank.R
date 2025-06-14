px_2PL_4_rank <- function(theta,a,d,xx=c(1,2,3,4)){ #a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)
  z <- exp(a*theta-d)
  prod(z[xx[-4]])/(sum(z)*sum(z[xx[-1]])*sum(z[xx[-1:-2]]))
}

grad_theta_2PL_4_rank <- function(theta,a,d,x=c(1,2,3,4)){
  z <- exp(a*theta-d)
  ai <- a[1]
  aj <- a[2]
  ak <- a[3]
  al <- a[4]
  # p.k.ijk <- p.k.ikj <- p.k.jki <- p.k.jik <- p.k.kij <- p.k.kji
  p.k <- z[3]/sum(z)
  # p.j.ijk <- p.j.ikj <- p.j.jki <- p.j.jik <- p.j.kij <- p.j.kji
  p.j <- z[2]/sum(z)
  # p.i.ijk <- p.i.ikj <- p.i.jki <- p.i.jik <- p.i.kij <- p.i.kji
  p.i <- z[1]/sum(z)
  p.l <- z[4]/sum(z)
  p.i.ijk <- p.i.ikj <- p.i.jik <- p.i.jki <- p.i.kij <- p.i.kji <- z[1]/sum(z[-4])
  p.j.ijk <- p.j.ikj <- p.j.jik <- p.j.jki <- p.j.kij <- p.j.kji <- z[2]/sum(z[-4])
  p.k.ijk <- p.k.ikj <- p.k.jik <- p.k.jki <- p.k.kij <- p.k.kji <- z[3]/sum(z[-4])
  p.i.ikl <- p.i.ilk <- p.i.kil <- p.i.kli <- p.i.lik <- p.i.lki <- z[1]/sum(z[-2])
  p.k.ikl <- p.k.ilk <- p.k.kil <- p.k.kli <- p.k.lik <- p.k.lki <- z[3]/sum(z[-2])
  p.l.ikl <- p.l.ilk <- p.l.kil <- p.l.kli <- p.l.lik <- p.l.lki <- z[4]/sum(z[-2])
  p.i.ijl <- p.i.ilj <- p.i.jil <- p.i.jli <- p.i.lij <- p.i.lji <- z[1]/sum(z[-3])
  p.j.ijl <- p.j.ilj <- p.j.jil <- p.j.jli <- p.j.lij <- p.j.lji <- z[2]/sum(z[-3])
  p.l.ijl <- p.l.ilj <- p.l.jil <- p.l.jli <- p.l.lij <- p.l.lji <- z[4]/sum(z[-3])
  p.j.jkl <- p.j.jlk <- p.j.klj <- p.j.kjl <- p.j.ljk <- p.j.lkj <- z[2]/sum(z[-1])
  p.k.jkl <- p.k.jlk <- p.k.klj <- p.k.kjl <- p.k.ljk <- p.k.lkj <- z[3]/sum(z[-1])
  p.l.jkl <- p.l.jlk <- p.l.klj <- p.l.kjl <- p.l.ljk <- p.l.lkj <- z[4]/sum(z[-1])

  p.k.kl <- p.k.lk <- z[3]/sum(z[-1:-2])
  p.l.kl <- p.l.lk <- z[4]/sum(z[-1:-2])
  p.i.ij <- p.i.ji <- z[1]/sum(z[-3:-4])
  p.j.ij <- p.j.ji <- z[2]/sum(z[-3:-4])
  p.i.ik <- p.i.ki <- z[1]/sum(z[c(-2,-4)])
  p.k.ik <- p.k.ki <- z[3]/sum(z[c(-2,-4)])
  p.i.il <- p.i.li <- z[1]/sum(z[-2:-3])
  p.l.il <- p.l.li <- z[4]/sum(z[-2:-3])
  p.j.jk <- p.j.kj <- z[2]/sum(z[c(-1,-4)])
  p.k.jk <- p.k.kj <- z[3]/sum(z[c(-1,-4)])
  p.j.jl <- p.j.lj <- z[2]/sum(z[c(-1,-3)])
  p.l.jl <- p.l.lj <- z[4]/sum(z[c(-1,-3)])
  if(identical(x,c(1,2,3,4))){
    g.i <- ai*p.i*(1-p.i)*p.j.jkl*p.k.kl
    g.j <- aj*p.i*p.j.jkl*p.k.kl*(1-p.j-p.j.jkl)
    g.k <- ak*p.i*p.j.jkl*p.k.kl*(p.l.kl-p.k.jkl-p.k)
    g.l <- -al*p.i*p.j.jkl*p.k.kl*(p.l+p.l.jkl+p.l.kl)
  }else if(identical(x,c(1,2,4,3))){
    g.i <- ai*p.i*(1-p.i)*p.j.jlk*p.l.lk
    g.j <- aj*p.i*p.j.jlk*p.l.lk*(1-p.j-p.j.jlk)
    g.l <- al*p.i*p.j.jlk*p.l.lk*(p.k.lk-p.l.jlk-p.l)
    g.k <- -ak*p.i*p.j.jlk*p.l.lk*(p.k+p.k.jlk+p.k.lk)
  }else if(identical(x,c(1,3,2,4))){
    g.i <- ai*p.i*(1-p.i)*p.k.kjl*p.j.jl
    g.k <- ak*p.i*p.k.kjl*p.j.jl*(1-p.k-p.k.kjl)
    g.j <- aj*p.i*p.k.kjl*p.j.jl*(p.l.jl-p.j.kjl-p.j)
    g.l <- -al*p.i*p.k.kjl*p.j.jl*(p.l+p.l.kjl+p.l.jl)
  }else if(identical(x,c(1,3,4,2))){
    g.i <- ai*p.i*(1-p.i)*p.k.klj*p.l.lj
    g.k <- ak*p.i*p.k.klj*p.l.lj*(1-p.k-p.k.klj)
    g.l <- al*p.i*p.k.klj*p.l.lj*(p.j.lj-p.l.klj-p.l)
    g.j <- -aj*p.i*p.k.klj*p.l.lj*(p.j+p.j.klj+p.j.lj)
  }else if(identical(x,c(1,4,2,3))){
    g.i <- ai*p.i*(1-p.i)*p.l.ljk*p.j.jk
    g.l <- al*p.i*p.l.ljk*p.j.jk*(1-p.l-p.l.ljk)
    g.j <- aj*p.i*p.l.ljk*p.j.jk*(p.k.jk-p.j.ljk-p.j)
    g.k <- -ak*p.i*p.l.ljk*p.j.jk*(p.k+p.k.ljk+p.k.jk)
  }else if(identical(x,c(1,4,3,2))){
    g.i <- ai*p.i*(1-p.i)*p.l.lkj*p.k.kj
    g.l <- al*p.i*p.l.lkj*p.k.kj*(1-p.l-p.l.lkj)
    g.k <- ak*p.i*p.l.lkj*p.k.kj*(p.j.kj-p.k.lkj-p.k)
    g.j <- -aj*p.i*p.l.lkj*p.k.kj*(p.j+p.j.lkj+p.j.kj)
  }else if(identical(x,c(2,1,3,4))){
    g.j <- aj*p.j*(1-p.j)*p.i.ikl*p.k.kl
    g.i <- ai*p.j*p.i.ikl*p.k.kl*(1-p.i-p.i.ikl)
    g.k <- ak*p.j*p.i.ikl*p.k.kl*(p.l.kl-p.k.ikl-p.k)
    g.l <- -al*p.j*p.i.ikl*p.k.kl*(p.l+p.l.ikl+p.l.kl)
  }else if(identical(x,c(2,1,4,3))){
    g.j <- aj*p.j*(1-p.j)*p.i.ilk*p.l.lk
    g.i <- ai*p.j*p.i.ilk*p.l.lk*(1-p.i-p.i.ilk)
    g.l <- al*p.j*p.i.ilk*p.l.lk*(p.k.lk-p.l.ilk-p.l)
    g.k <- -ak*p.j*p.i.ilk*p.l.lk*(p.k+p.k.ilk+p.k.lk)
  }else if(identical(x,c(2,3,1,4))){
    g.j <- aj*p.j*(1-p.j)*p.k.kil*p.i.il
    g.k <- ak*p.j*p.k.kil*p.i.il*(1-p.k-p.k.kil)
    g.i <- ai*p.j*p.k.kil*p.i.il*(p.l.il-p.i.kil-p.i)
    g.l <- -al*p.j*p.k.kil*p.i.il*(p.l+p.l.kil+p.l.il)
  }else if(identical(x,c(2,3,4,1))){
    g.j <- aj*p.j*(1-p.j)*p.k.kli*p.l.li
    g.k <- ak*p.j*p.k.kli*p.l.li*(1-p.k-p.k.kli)
    g.l <- al*p.j*p.k.kli*p.l.li*(p.i.li-p.l.kli-p.l)
    g.i <- -ai*p.j*p.k.kli*p.l.li*(p.i+p.i.kli+p.i.li)
  }else if(identical(x,c(2,4,1,3))){
    g.j <- aj*p.j*(1-p.j)*p.l.lik*p.i.ik
    g.l <- al*p.j*p.l.lik*p.i.ik*(1-p.l-p.l.lik)
    g.i <- ai*p.j*p.l.lik*p.i.ik*(p.k.ik-p.i.lik-p.i)
    g.k <- -ak*p.j*p.l.lik*p.i.ik*(p.k+p.k.lik+p.k.ik)
  }else if(identical(x,c(2,4,3,1))){
    g.j <- aj*p.j*(1-p.j)*p.l.lki*p.k.ki
    g.l <- al*p.j*p.l.lki*p.k.ki*(1-p.l-p.l.lki)
    g.k <- ak*p.j*p.l.lki*p.k.ki*(p.i.ki-p.k.lki-p.k)
    g.i <- -ai*p.j*p.l.lki*p.k.ki*(p.i+p.i.lki+p.i.ki)
  }else if(identical(x,c(3,1,2,4))){
    g.k <- ak*p.k*(1-p.k)*p.i.ijl*p.j.jl
    g.i <- ai*p.k*p.i.ijl*p.j.jl*(1-p.i-p.i.ijl)
    g.j <- aj*p.k*p.i.ijl*p.j.jl*(p.l.jl-p.j.ijl-p.j)
    g.l <- -al*p.k*p.i.ijl*p.j.jl*(p.l+p.l.ijl+p.l.jl)
  }else if(identical(x,c(3,1,4,2))){
    g.k <- ak*p.k*(1-p.k)*p.i.ilj*p.l.lj
    g.i <- ai*p.k*p.i.ilj*p.l.lj*(1-p.i-p.i.ilj)
    g.l <- al*p.k*p.i.ilj*p.l.lj*(p.j.lj-p.l.ilj-p.l)
    g.j <- -aj*p.k*p.i.ilj*p.l.lj*(p.j+p.j.ilj+p.j.lj)
  }else if(identical(x,c(3,2,1,4))){
    g.k <- ak*p.k*(1-p.k)*p.j.jil*p.i.il
    g.j <- aj*p.k*p.j.jil*p.i.il*(1-p.j-p.j.jil)
    g.i <- ai*p.k*p.j.jil*p.i.il*(p.l.il-p.i.jil-p.i)
    g.l <- -al*p.k*p.j.jil*p.i.il*(p.l+p.l.jil+p.l.il)
  }else if(identical(x,c(3,2,4,1))){
    g.k <- ak*p.k*(1-p.k)*p.j.jli*p.l.li
    g.j <- aj*p.k*p.j.jli*p.l.li*(1-p.j-p.j.jli)
    g.l <- al*p.k*p.j.jli*p.l.li*(p.i.li-p.l.jli-p.l)
    g.i <- -ai*p.k*p.j.jli*p.l.li*(p.i+p.i.jli+p.i.li)
  }else if(identical(x,c(3,4,1,2))){
    g.k <- ak*p.k*(1-p.k)*p.l.lij*p.i.ij
    g.l <- al*p.k*p.l.lij*p.i.ij*(1-p.l-p.l.lij)
    g.i <- ai*p.k*p.l.lij*p.i.ij*(p.j.ij-p.i.lij-p.i)
    g.j <- -aj*p.k*p.l.lij*p.i.ij*(p.j+p.j.lij+p.j.ij)
  }else if(identical(x,c(3,4,2,1))){
    g.k <- ak*p.k*(1-p.k)*p.l.lji*p.j.ji
    g.l <- al*p.k*p.l.lji*p.j.ji*(1-p.l-p.l.lji)
    g.j <- aj*p.k*p.l.lji*p.j.ji*(p.i.ji-p.j.lji-p.j)
    g.i <- -ai*p.k*p.l.lji*p.j.ji*(p.i+p.i.lji+p.i.ji)
  }else if(identical(x,c(4,1,2,3))){
    g.l <- al*p.l*(1-p.l)*p.i.ijk*p.j.jk
    g.i <- ai*p.l*p.i.ijk*p.j.jk*(1-p.i-p.i.ijk)
    g.j <- aj*p.l*p.i.ijk*p.j.jk*(p.k.jk-p.j.ijk-p.j)
    g.k <- -ak*p.l*p.i.ijk*p.j.jk*(p.k+p.k.ijk+p.k.jk)
  }else if(identical(x,c(4,1,3,2))){
    g.l <- al*p.l*(1-p.l)*p.i.ikj*p.k.kj
    g.i <- ai*p.l*p.i.ikj*p.k.kj*(1-p.i-p.i.ikj)
    g.k <- ak*p.l*p.i.ikj*p.k.kj*(p.j.kj-p.k.ikj-p.k)
    g.j <- -aj*p.l*p.i.ikj*p.k.kj*(p.j+p.j.ikj+p.j.kj)
  }else if(identical(x,c(4,2,1,3))){
    g.l <- al*p.l*(1-p.l)*p.j.jik*p.i.ik
    g.j <- aj*p.l*p.j.jik*p.i.ik*(1-p.j-p.j.jik)
    g.i <- ai*p.l*p.j.jik*p.i.ik*(p.k.ik-p.i.jik-p.i)
    g.k <- -ak*p.l*p.j.jik*p.i.ik*(p.k+p.k.jik+p.k.ik)
  }else if(identical(x,c(4,2,3,1))){
    g.l <- al*p.l*(1-p.l)*p.j.jki*p.k.ki
    g.j <- aj*p.l*p.j.jki*p.k.ki*(1-p.j-p.j.jki)
    g.k <- ak*p.l*p.j.jki*p.k.ki*(p.i.ki-p.k.jki-p.k)
    g.i <- -ai*p.l*p.j.jki*p.k.ki*(p.i+p.i.jki+p.i.ki)
  }else if(identical(x,c(4,3,1,2))){
    g.l <- al*p.l*(1-p.l)*p.k.kij*p.i.ij
    g.k <- ak*p.l*p.k.kij*p.i.ij*(1-p.k-p.k.kij)
    g.i <- ai*p.l*p.k.kij*p.i.ij*(p.j.ij-p.i.kij-p.i)
    g.j <- -aj*p.l*p.k.kij*p.i.ij*(p.j+p.j.kij+p.j.ij)
  }else if(identical(x,c(4,3,2,1))){
    g.l <- al*p.l*(1-p.l)*p.k.kji*p.j.ji
    g.k <- ak*p.l*p.k.kji*p.j.ji*(1-p.k-p.k.kji)
    g.j <- aj*p.l*p.k.kji*p.j.ji*(p.i.ji-p.j.kji-p.j)
    g.i <- -ai*p.l*p.k.kji*p.j.ji*(p.i+p.i.kji+p.i.ji)
  }

  c(g.i,g.j,g.k,g.l)
}

hes_theta_2PL_4_rank <- function(theta,a,d,x=c(1,2,3,4)){
  h <- matrix(NA,4,4)
  z <- exp(a*theta-d)
  ai <- a[1]
  aj <- a[2]
  ak <- a[3]
  al <- a[4]
  # p.k.ijk <- p.k.ikj <- p.k.jki <- p.k.jik <- p.k.kij <- p.k.kji
  p.k <- z[3]/sum(z)
  p.l <- z[4]/sum(z)
  # p.j.ijk <- p.j.ikj <- p.j.jki <- p.j.jik <- p.j.kij <- p.j.kji
  p.j <- z[2]/sum(z)
  # p.i.ijk <- p.i.ikj <- p.i.jki <- p.i.jik <- p.i.kij <- p.i.kji
  p.i <- z[1]/sum(z)
  p.i.ijk <- p.i.ikj <- p.i.jik <- p.i.jki <- p.i.kij <- p.i.kji <- z[1]/sum(z[-4])
  p.j.ijk <- p.j.ikj <- p.j.jik <- p.j.jki <- p.j.kij <- p.j.kji <- z[2]/sum(z[-4])
  p.k.ijk <- p.k.ikj <- p.k.jik <- p.k.jki <- p.k.kij <- p.k.kji <- z[3]/sum(z[-4])
  p.i.ikl <- p.i.ilk <- p.i.kil <- p.i.kli <- p.i.lik <- p.i.lki <- z[1]/sum(z[-2])
  p.k.ikl <- p.k.ilk <- p.k.kil <- p.k.kli <- p.k.lik <- p.k.lki <- z[3]/sum(z[-2])
  p.l.ikl <- p.l.ilk <- p.l.kil <- p.l.kli <- p.l.lik <- p.l.lki <- z[4]/sum(z[-2])
  p.i.ijl <- p.i.ilj <- p.i.jil <- p.i.jli <- p.i.lij <- p.i.lji <- z[1]/sum(z[-3])
  p.j.ijl <- p.j.ilj <- p.j.jil <- p.j.jli <- p.j.lij <- p.j.lji <- z[2]/sum(z[-3])
  p.l.ijl <- p.l.ilj <- p.l.jil <- p.l.jli <- p.l.lij <- p.l.lji <- z[4]/sum(z[-3])
  p.j.jkl <- p.j.jlk <- p.j.klj <- p.j.kjl <- p.j.ljk <- p.j.lkj <- z[2]/sum(z[-1])
  p.k.jkl <- p.k.jlk <- p.k.klj <- p.k.kjl <- p.k.ljk <- p.k.lkj <- z[3]/sum(z[-1])
  p.l.jkl <- p.l.jlk <- p.l.klj <- p.l.kjl <- p.l.ljk <- p.l.lkj <- z[4]/sum(z[-1])

  p.k.kl <- p.k.lk <- z[3]/sum(z[-1:-2])
  p.l.kl <- p.l.lk <- z[4]/sum(z[-1:-2])
  p.i.ij <- p.i.ji <- z[1]/sum(z[-3:-4])
  p.j.ij <- p.j.ji <- z[2]/sum(z[-3:-4])
  p.i.ik <- p.i.ki <- z[1]/sum(z[c(-2,-4)])
  p.k.ik <- p.k.ki <- z[3]/sum(z[c(-2,-4)])
  p.i.il <- p.i.li <- z[1]/sum(z[-2:-3])
  p.l.il <- p.l.li <- z[4]/sum(z[-2:-3])
  p.j.jk <- p.j.kj <- z[2]/sum(z[c(-1,-4)])
  p.k.jk <- p.k.kj <- z[3]/sum(z[c(-1,-4)])
  p.j.jl <- p.j.lj <- z[2]/sum(z[c(-1,-3)])
  p.l.jl <- p.l.lj <- z[4]/sum(z[c(-1,-3)])


    if(identical(x,c(1,2,3,4))){
      ii <- ai^2*p.i*p.j.jkl*(1-2*p.i)*p.k.kl*(1-p.i)
      ij <- ji <- ai*aj*p.i*p.j.jkl*(1-2*p.j)*p.k.kl*(1-p.i)
      ik <- ki <- ai*ak*p.i*p.j.jkl*p.k.kl*(2*p.k*p.i-2*p.k+(1-p.i)*p.l.kl)
      il <- li <- ai*al*p.i*p.j.jkl*p.k.kl*(2*p.l*p.i-2*p.l-(1-p.i)*p.l.kl)
      jj <- aj^2*p.i*p.j.jkl*p.k.kl*(1-3*p.j+2*p.j^2+2*p.j*p.j.jkl+2*p.j.jkl^2)
      jk <- kj <- aj*ak*p.i*p.j.jkl*p.k.kl*((p.l.kl-p.k-p.k.jkl)*(1-p.j-p.j.jkl)+p.j*p.k+p.j.jkl*p.k.jkl)
      jl <- lj <- aj*al*p.i*p.j.jkl*p.k.kl*((-p.l.kl-p.l-p.l.jkl)*(1-p.j-p.j.jkl)+p.j*p.l+p.j.jkl*p.l.jkl)
      kk <- ak^2*p.i*p.j.jkl*p.k.kl*((p.l.kl-p.k-p.k.jkl)^2-p.l.kl*p.k.kl-p.k.jkl*(1-p.k.jkl)-p.k*(1-p.k))
      kl <- lk <- ak*al*p.i*p.j.jkl*p.k.kl*((p.k+p.k.jkl-p.l.kl)*(p.l+p.l.jkl+p.l.kl)+p.k.jkl*p.l.jkl+p.l*p.k+p.k.kl*p.l.kl)
      ll <- al^2*p.i*p.j.jkl*p.k.kl*((p.l+p.l.jkl+p.l.kl)^2-p.l.jkl*(1-p.l.jkl)-p.l*(1-p.l)-p.l.kl*(1-p.l.kl))
    }else if(identical(x,c(1,2,4,3))){
      ii <- ai^2*p.i*p.j.jlk*(1-2*p.i)*p.l.lk*(1-p.i)
      ij <- ji <- ai*aj*p.i*p.j.jlk*(1-2*p.j)*p.l.lk*(1-p.i)
      il <- li <- ai*al*p.i*p.j.jlk*p.l.lk*(2*p.l*p.i-2*p.l+(1-p.i)*p.k.lk)
      ik <- ki <- ai*ak*p.i*p.j.jlk*p.l.lk*(2*p.k*p.i-2*p.k-(1-p.i)*p.k.lk)
      jj <- aj^2*p.i*p.j.jlk*p.l.lk*(1-3*p.j+2*p.j^2+2*p.j*p.j.jlk+2*p.j.jlk^2)
      jl <- lj <- aj*al*p.i*p.j.jlk*p.l.lk*((p.k.lk-p.l-p.l.jlk)*(1-p.j-p.j.jlk)+p.j*p.l+p.j.jlk*p.l.jlk)
      jk <- kj <- aj*ak*p.i*p.j.jlk*p.l.lk*((-p.k.lk-p.k-p.k.jlk)*(1-p.j-p.j.jlk)+p.j*p.k+p.j.jlk*p.k.jlk)
      ll <- al^2*p.i*p.j.jlk*p.l.lk*((p.k.lk-p.l-p.l.jlk)^2-p.k.lk*p.l.lk-p.l.jlk*(1-p.l.jlk)-p.l*(1-p.l))
      lk <- kl <- al*ak*p.i*p.j.jlk*p.l.lk*((p.l+p.l.jlk-p.k.lk)*(p.k+p.k.jlk+p.k.lk)+p.l.jlk*p.k.jlk+p.k*p.l+p.l.lk*p.k.lk)
      kk <- ak^2*p.i*p.j.jlk*p.l.lk*((p.k+p.k.jlk+p.k.lk)^2-p.k.jlk*(1-p.k.jlk)-p.k*(1-p.k)-p.k.lk*(1-p.k.lk))
    }else if(identical(x,c(1,3,2,4))){
      ii <- ai^2*p.i*p.k.kjl*(1-2*p.i)*p.j.jl*(1-p.i)
      ik <- ki <- ai*ak*p.i*p.k.kjl*(1-2*p.k)*p.j.jl*(1-p.i)
      ij <- ji <- ai*aj*p.i*p.k.kjl*p.j.jl*(2*p.j*p.i-2*p.j+(1-p.i)*p.l.jl)
      il <- li <- ai*al*p.i*p.k.kjl*p.j.jl*(2*p.l*p.i-2*p.l-(1-p.i)*p.l.jl)
      kk <- ak^2*p.i*p.k.kjl*p.j.jl*(1-3*p.k+2*p.k^2+2*p.k*p.k.kjl+2*p.k.kjl^2)
      kj <- jk <- ak*aj*p.i*p.k.kjl*p.j.jl*((p.l.jl-p.j-p.j.kjl)*(1-p.k-p.k.kjl)+p.k*p.j+p.k.kjl*p.j.kjl)
      kl <- lk <- ak*al*p.i*p.k.kjl*p.j.jl*((-p.l.jl-p.l-p.l.kjl)*(1-p.k-p.k.kjl)+p.k*p.l+p.k.kjl*p.l.kjl)
      jj <- aj^2*p.i*p.k.kjl*p.j.jl*((p.l.jl-p.j-p.j.kjl)^2-p.l.jl*p.j.jl-p.j.kjl*(1-p.j.kjl)-p.j*(1-p.j))
      jl <- lj <- aj*al*p.i*p.k.kjl*p.j.jl*((p.j+p.j.kjl-p.l.jl)*(p.l+p.l.kjl+p.l.jl)+p.j.kjl*p.l.kjl+p.l*p.j+p.j.jl*p.l.jl)
      ll <- al^2*p.i*p.k.kjl*p.j.jl*((p.l+p.l.kjl+p.l.jl)^2-p.l.kjl*(1-p.l.kjl)-p.l*(1-p.l)-p.l.jl*(1-p.l.jl))
    }else if(identical(x,c(1,3,4,2))){
      ii <- ai^2*p.i*p.k.klj*(1-2*p.i)*p.l.lj*(1-p.i)
      ik <- ki <- ai*ak*p.i*p.k.klj*(1-2*p.k)*p.l.lj*(1-p.i)
      il <- li <- ai*al*p.i*p.k.klj*p.l.lj*(2*p.l*p.i-2*p.l+(1-p.i)*p.j.lj)
      ij <- ji <- ai*aj*p.i*p.k.klj*p.l.lj*(2*p.j*p.i-2*p.j-(1-p.i)*p.j.lj)
      kk <- ak^2*p.i*p.k.klj*p.l.lj*(1-3*p.k+2*p.k^2+2*p.k*p.k.klj+2*p.k.klj^2)
      kl <- lk <- ak*al*p.i*p.k.klj*p.l.lj*((p.j.lj-p.l-p.l.klj)*(1-p.k-p.k.klj)+p.k*p.l+p.k.klj*p.l.klj)
      kj <- jk <- ak*aj*p.i*p.k.klj*p.l.lj*((-p.j.lj-p.j-p.j.klj)*(1-p.k-p.k.klj)+p.k*p.j+p.k.klj*p.j.klj)
      ll <- al^2*p.i*p.k.klj*p.l.lj*((p.j.lj-p.l-p.l.klj)^2-p.j.lj*p.l.lj-p.l.klj*(1-p.l.klj)-p.l*(1-p.l))
      lj <- jl <- al*aj*p.i*p.k.klj*p.l.lj*((p.l+p.l.klj-p.j.lj)*(p.j+p.j.klj+p.j.lj)+p.l.klj*p.j.klj+p.j*p.l+p.l.lj*p.j.lj)
      jj <- aj^2*p.i*p.k.klj*p.l.lj*((p.j+p.j.klj+p.j.lj)^2-p.j.klj*(1-p.j.klj)-p.j*(1-p.j)-p.j.lj*(1-p.j.lj))
    }else if(identical(x,c(1,4,2,3))){
      ii <- ai^2*p.i*p.l.ljk*(1-2*p.i)*p.j.jk*(1-p.i)
      il <- li <- ai*al*p.i*p.l.ljk*(1-2*p.l)*p.j.jk*(1-p.i)
      ij <- ji <- ai*aj*p.i*p.l.ljk*p.j.jk*(2*p.j*p.i-2*p.j+(1-p.i)*p.k.jk)
      ik <- ki <- ai*ak*p.i*p.l.ljk*p.j.jk*(2*p.k*p.i-2*p.k-(1-p.i)*p.k.jk)
      ll <- al^2*p.i*p.l.ljk*p.j.jk*(1-3*p.l+2*p.l^2+2*p.l*p.l.ljk+2*p.l.ljk^2)
      lj <- jl <- al*aj*p.i*p.l.ljk*p.j.jk*((p.k.jk-p.j-p.j.ljk)*(1-p.l-p.l.ljk)+p.l*p.j+p.l.ljk*p.j.ljk)
      lk <- kl <- al*ak*p.i*p.l.ljk*p.j.jk*((-p.k.jk-p.k-p.k.ljk)*(1-p.l-p.l.ljk)+p.l*p.k+p.l.ljk*p.k.ljk)
      jj <- aj^2*p.i*p.l.ljk*p.j.jk*((p.k.jk-p.j-p.j.ljk)^2-p.k.jk*p.j.jk-p.j.ljk*(1-p.j.ljk)-p.j*(1-p.j))
      jk <- kj <- aj*ak*p.i*p.l.ljk*p.j.jk*((p.j+p.j.ljk-p.k.jk)*(p.k+p.k.ljk+p.k.jk)+p.j.ljk*p.k.ljk+p.k*p.j+p.j.jk*p.k.jk)
      kk <- ak^2*p.i*p.l.ljk*p.j.jk*((p.k+p.k.ljk+p.k.jk)^2-p.k.ljk*(1-p.k.ljk)-p.k*(1-p.k)-p.k.jk*(1-p.k.jk))
    }else if(identical(x,c(1,4,3,2))){
      ii <- ai^2*p.i*p.l.lkj*(1-2*p.i)*p.k.kj*(1-p.i)
      il <- li <- ai*al*p.i*p.l.lkj*(1-2*p.l)*p.k.kj*(1-p.i)
      ik <- ki <- ai*ak*p.i*p.l.lkj*p.k.kj*(2*p.k*p.i-2*p.k+(1-p.i)*p.j.kj)
      ij <- ji <- ai*aj*p.i*p.l.lkj*p.k.kj*(2*p.j*p.i-2*p.j-(1-p.i)*p.j.kj)
      ll <- al^2*p.i*p.l.lkj*p.k.kj*(1-3*p.l+2*p.l^2+2*p.l*p.l.lkj+2*p.l.lkj^2)
      lk <- kl <- al*ak*p.i*p.l.lkj*p.k.kj*((p.j.kj-p.k-p.k.lkj)*(1-p.l-p.l.lkj)+p.l*p.k+p.l.lkj*p.k.lkj)
      lj <- jl <- al*aj*p.i*p.l.lkj*p.k.kj*((-p.j.kj-p.j-p.j.lkj)*(1-p.l-p.l.lkj)+p.l*p.j+p.l.lkj*p.j.lkj)
      kk <- ak^2*p.i*p.l.lkj*p.k.kj*((p.j.kj-p.k-p.k.lkj)^2-p.j.kj*p.k.kj-p.k.lkj*(1-p.k.lkj)-p.k*(1-p.k))
      kj <- jk <- ak*aj*p.i*p.l.lkj*p.k.kj*((p.k+p.k.lkj-p.j.kj)*(p.j+p.j.lkj+p.j.kj)+p.k.lkj*p.j.lkj+p.j*p.k+p.k.kj*p.j.kj)
      jj <- aj^2*p.i*p.l.lkj*p.k.kj*((p.j+p.j.lkj+p.j.kj)^2-p.j.lkj*(1-p.j.lkj)-p.j*(1-p.j)-p.j.kj*(1-p.j.kj))
    }else if(identical(x,c(2,1,3,4))){
      jj <- aj^2*p.j*p.i.ikl*(1-2*p.j)*p.k.kl*(1-p.j)
      ji <- ij <- aj*ai*p.j*p.i.ikl*(1-2*p.i)*p.k.kl*(1-p.j)
      jk <- kj <- aj*ak*p.j*p.i.ikl*p.k.kl*(2*p.k*p.j-2*p.k+(1-p.j)*p.l.kl)
      jl <- lj <- aj*al*p.j*p.i.ikl*p.k.kl*(2*p.l*p.j-2*p.l-(1-p.j)*p.l.kl)
      ii <- ai^2*p.j*p.i.ikl*p.k.kl*(1-3*p.i+2*p.i^2+2*p.i*p.i.ikl+2*p.i.ikl^2)
      ik <- ki <- ai*ak*p.j*p.i.ikl*p.k.kl*((p.l.kl-p.k-p.k.ikl)*(1-p.i-p.i.ikl)+p.i*p.k+p.i.ikl*p.k.ikl)
      il <- li <- ai*al*p.j*p.i.ikl*p.k.kl*((-p.l.kl-p.l-p.l.ikl)*(1-p.i-p.i.ikl)+p.i*p.l+p.i.ikl*p.l.ikl)
      kk <- ak^2*p.j*p.i.ikl*p.k.kl*((p.l.kl-p.k-p.k.ikl)^2-p.l.kl*p.k.kl-p.k.ikl*(1-p.k.ikl)-p.k*(1-p.k))
      kl <- lk <- ak*al*p.j*p.i.ikl*p.k.kl*((p.k+p.k.ikl-p.l.kl)*(p.l+p.l.ikl+p.l.kl)+p.k.ikl*p.l.ikl+p.l*p.k+p.k.kl*p.l.kl)
      ll <- al^2*p.j*p.i.ikl*p.k.kl*((p.l+p.l.ikl+p.l.kl)^2-p.l.ikl*(1-p.l.ikl)-p.l*(1-p.l)-p.l.kl*(1-p.l.kl))
    }else if(identical(x,c(2,1,4,3))){
      jj <- aj^2*p.j*p.i.ilk*(1-2*p.j)*p.l.lk*(1-p.j)
      ji <- ij <- aj*ai*p.j*p.i.ilk*(1-2*p.i)*p.l.lk*(1-p.j)
      jl <- lj <- aj*al*p.j*p.i.ilk*p.l.lk*(2*p.l*p.j-2*p.l+(1-p.j)*p.k.lk)
      jk <- kj <- aj*ak*p.j*p.i.ilk*p.l.lk*(2*p.k*p.j-2*p.k-(1-p.j)*p.k.lk)
      ii <- ai^2*p.j*p.i.ilk*p.l.lk*(1-3*p.i+2*p.i^2+2*p.i*p.i.ilk+2*p.i.ilk^2)
      il <- li <- ai*al*p.j*p.i.ilk*p.l.lk*((p.k.lk-p.l-p.l.ilk)*(1-p.i-p.i.ilk)+p.i*p.l+p.i.ilk*p.l.ilk)
      ik <- ki <- ai*ak*p.j*p.i.ilk*p.l.lk*((-p.k.lk-p.k-p.k.ilk)*(1-p.i-p.i.ilk)+p.i*p.k+p.i.ilk*p.k.ilk)
      ll <- al^2*p.j*p.i.ilk*p.l.lk*((p.k.lk-p.l-p.l.ilk)^2-p.k.lk*p.l.lk-p.l.ilk*(1-p.l.ilk)-p.l*(1-p.l))
      lk <- kl <- al*ak*p.j*p.i.ilk*p.l.lk*((p.l+p.l.ilk-p.k.lk)*(p.k+p.k.ilk+p.k.lk)+p.l.ilk*p.k.ilk+p.k*p.l+p.l.lk*p.k.lk)
      kk <- ak^2*p.j*p.i.ilk*p.l.lk*((p.k+p.k.ilk+p.k.lk)^2-p.k.ilk*(1-p.k.ilk)-p.k*(1-p.k)-p.k.lk*(1-p.k.lk))
    }else if(identical(x,c(2,3,1,4))){
      jj <- aj^2*p.j*p.k.kil*(1-2*p.j)*p.i.il*(1-p.j)
      jk <- kj <- aj*ak*p.j*p.k.kil*(1-2*p.k)*p.i.il*(1-p.j)
      ji <- ij <- aj*ai*p.j*p.k.kil*p.i.il*(2*p.i*p.j-2*p.i+(1-p.j)*p.l.il)
      jl <- lj <- aj*al*p.j*p.k.kil*p.i.il*(2*p.l*p.j-2*p.l-(1-p.j)*p.l.il)
      kk <- ak^2*p.j*p.k.kil*p.i.il*(1-3*p.k+2*p.k^2+2*p.k*p.k.kil+2*p.k.kil^2)
      ki <- ik <- ak*ai*p.j*p.k.kil*p.i.il*((p.l.il-p.i-p.i.kil)*(1-p.k-p.k.kil)+p.k*p.i+p.k.kil*p.i.kil)
      kl <- lk <- ak*al*p.j*p.k.kil*p.i.il*((-p.l.il-p.l-p.l.kil)*(1-p.k-p.k.kil)+p.k*p.l+p.k.kil*p.l.kil)
      ii <- ai^2*p.j*p.k.kil*p.i.il*((p.l.il-p.i-p.i.kil)^2-p.l.il*p.i.il-p.i.kil*(1-p.i.kil)-p.i*(1-p.i))
      il <- li <- ai*al*p.j*p.k.kil*p.i.il*((p.i+p.i.kil-p.l.il)*(p.l+p.l.kil+p.l.il)+p.i.kil*p.l.kil+p.l*p.i+p.i.il*p.l.il)
      ll <- al^2*p.j*p.k.kil*p.i.il*((p.l+p.l.kil+p.l.il)^2-p.l.kil*(1-p.l.kil)-p.l*(1-p.l)-p.l.il*(1-p.l.il))
    }else if(identical(x,c(2,3,4,1))){
      jj <- aj^2*p.j*p.k.kli*(1-2*p.j)*p.l.li*(1-p.j)
      jk <- kj <- aj*ak*p.j*p.k.kli*(1-2*p.k)*p.l.li*(1-p.j)
      jl <- lj <- aj*al*p.j*p.k.kli*p.l.li*(2*p.l*p.j-2*p.l+(1-p.j)*p.i.li)
      ji <- ij <- aj*ai*p.j*p.k.kli*p.l.li*(2*p.i*p.j-2*p.i-(1-p.j)*p.i.li)
      kk <- ak^2*p.j*p.k.kli*p.l.li*(1-3*p.k+2*p.k^2+2*p.k*p.k.kli+2*p.k.kli^2)
      kl <- lk <- ak*al*p.j*p.k.kli*p.l.li*((p.i.li-p.l-p.l.kli)*(1-p.k-p.k.kli)+p.k*p.l+p.k.kli*p.l.kli)
      ki <- ik <- ak*ai*p.j*p.k.kli*p.l.li*((-p.i.li-p.i-p.i.kli)*(1-p.k-p.k.kli)+p.k*p.i+p.k.kli*p.i.kli)
      ll <- al^2*p.j*p.k.kli*p.l.li*((p.i.li-p.l-p.l.kli)^2-p.i.li*p.l.li-p.l.kli*(1-p.l.kli)-p.l*(1-p.l))
      li <- il <- al*ai*p.j*p.k.kli*p.l.li*((p.l+p.l.kli-p.i.li)*(p.i+p.i.kli+p.i.li)+p.l.kli*p.i.kli+p.i*p.l+p.l.li*p.i.li)
      ii <- ai^2*p.j*p.k.kli*p.l.li*((p.i+p.i.kli+p.i.li)^2-p.i.kli*(1-p.i.kli)-p.i*(1-p.i)-p.i.li*(1-p.i.li))
    }else if(identical(x,c(2,4,1,3))){
      jj <- aj^2*p.j*p.l.lik*(1-2*p.j)*p.i.ik*(1-p.j)
      jl <- lj <- aj*al*p.j*p.l.lik*(1-2*p.l)*p.i.ik*(1-p.j)
      ji <- ij <- aj*ai*p.j*p.l.lik*p.i.ik*(2*p.i*p.j-2*p.i+(1-p.j)*p.k.ik)
      jk <- kj <- aj*ak*p.j*p.l.lik*p.i.ik*(2*p.k*p.j-2*p.k-(1-p.j)*p.k.ik)
      ll <- al^2*p.j*p.l.lik*p.i.ik*(1-3*p.l+2*p.l^2+2*p.l*p.l.lik+2*p.l.lik^2)
      li <- il <- al*ai*p.j*p.l.lik*p.i.ik*((p.k.ik-p.i-p.i.lik)*(1-p.l-p.l.lik)+p.l*p.i+p.l.lik*p.i.lik)
      lk <- kl <- al*ak*p.j*p.l.lik*p.i.ik*((-p.k.ik-p.k-p.k.lik)*(1-p.l-p.l.lik)+p.l*p.k+p.l.lik*p.k.lik)
      ii <- ai^2*p.j*p.l.lik*p.i.ik*((p.k.ik-p.i-p.i.lik)^2-p.k.ik*p.i.ik-p.i.lik*(1-p.i.lik)-p.i*(1-p.i))
      ik <- ki <- ai*ak*p.j*p.l.lik*p.i.ik*((p.i+p.i.lik-p.k.ik)*(p.k+p.k.lik+p.k.ik)+p.i.lik*p.k.lik+p.k*p.i+p.i.ik*p.k.ik)
      kk <- ak^2*p.j*p.l.lik*p.i.ik*((p.k+p.k.lik+p.k.ik)^2-p.k.lik*(1-p.k.lik)-p.k*(1-p.k)-p.k.ik*(1-p.k.ik))
    }else if(identical(x,c(2,4,3,1))){
      jj <- aj^2*p.j*p.l.lki*(1-2*p.j)*p.k.ki*(1-p.j)
      jl <- lj <- aj*al*p.j*p.l.lki*(1-2*p.l)*p.k.ki*(1-p.j)
      jk <- kj <- aj*ak*p.j*p.l.lki*p.k.ki*(2*p.k*p.j-2*p.k+(1-p.j)*p.i.ki)
      ji <- ij <- aj*ai*p.j*p.l.lki*p.k.ki*(2*p.i*p.j-2*p.i-(1-p.j)*p.i.ki)
      ll <- al^2*p.j*p.l.lki*p.k.ki*(1-3*p.l+2*p.l^2+2*p.l*p.l.lki+2*p.l.lki^2)
      lk <- kl <- al*ak*p.j*p.l.lki*p.k.ki*((p.i.ki-p.k-p.k.lki)*(1-p.l-p.l.lki)+p.l*p.k+p.l.lki*p.k.lki)
      li <- il <- al*ai*p.j*p.l.lki*p.k.ki*((-p.i.ki-p.i-p.i.lki)*(1-p.l-p.l.lki)+p.l*p.i+p.l.lki*p.i.lki)
      kk <- ak^2*p.j*p.l.lki*p.k.ki*((p.i.ki-p.k-p.k.lki)^2-p.i.ki*p.k.ki-p.k.lki*(1-p.k.lki)-p.k*(1-p.k))
      ki <- ik <- ak*ai*p.j*p.l.lki*p.k.ki*((p.k+p.k.lki-p.i.ki)*(p.i+p.i.lki+p.i.ki)+p.k.lki*p.i.lki+p.i*p.k+p.k.ki*p.i.ki)
      ii <- ai^2*p.j*p.l.lki*p.k.ki*((p.i+p.i.lki+p.i.ki)^2-p.i.lki*(1-p.i.lki)-p.i*(1-p.i)-p.i.ki*(1-p.i.ki))
    }else if(identical(x,c(3,1,2,4))){
      kk <- ak^2*p.k*p.i.ijl*(1-2*p.k)*p.j.jl*(1-p.k)
      ki <- ik <- ak*ai*p.k*p.i.ijl*(1-2*p.i)*p.j.jl*(1-p.k)
      kj <- jk <- ak*aj*p.k*p.i.ijl*p.j.jl*(2*p.j*p.k-2*p.j+(1-p.k)*p.l.jl)
      kl <- lk <- ak*al*p.k*p.i.ijl*p.j.jl*(2*p.l*p.k-2*p.l-(1-p.k)*p.l.jl)
      ii <- ai^2*p.k*p.i.ijl*p.j.jl*(1-3*p.i+2*p.i^2+2*p.i*p.i.ijl+2*p.i.ijl^2)
      ij <- ji <- ai*aj*p.k*p.i.ijl*p.j.jl*((p.l.jl-p.j-p.j.ijl)*(1-p.i-p.i.ijl)+p.i*p.j+p.i.ijl*p.j.ijl)
      il <- li <- ai*al*p.k*p.i.ijl*p.j.jl*((-p.l.jl-p.l-p.l.ijl)*(1-p.i-p.i.ijl)+p.i*p.l+p.i.ijl*p.l.ijl)
      jj <- aj^2*p.k*p.i.ijl*p.j.jl*((p.l.jl-p.j-p.j.ijl)^2-p.l.jl*p.j.jl-p.j.ijl*(1-p.j.ijl)-p.j*(1-p.j))
      jl <- lj <- aj*al*p.k*p.i.ijl*p.j.jl*((p.j+p.j.ijl-p.l.jl)*(p.l+p.l.ijl+p.l.jl)+p.j.ijl*p.l.ijl+p.l*p.j+p.j.jl*p.l.jl)
      ll <- al^2*p.k*p.i.ijl*p.j.jl*((p.l+p.l.ijl+p.l.jl)^2-p.l.ijl*(1-p.l.ijl)-p.l*(1-p.l)-p.l.jl*(1-p.l.jl))
    }else if(identical(x,c(3,1,4,2))){
      kk <- ak^2*p.k*p.i.ilj*(1-2*p.k)*p.l.lj*(1-p.k)
      ki <- ik <- ak*ai*p.k*p.i.ilj*(1-2*p.i)*p.l.lj*(1-p.k)
      kl <- lk <- ak*al*p.k*p.i.ilj*p.l.lj*(2*p.l*p.k-2*p.l+(1-p.k)*p.j.lj)
      kj <- jk <- ak*aj*p.k*p.i.ilj*p.l.lj*(2*p.j*p.k-2*p.j-(1-p.k)*p.j.lj)
      ii <- ai^2*p.k*p.i.ilj*p.l.lj*(1-3*p.i+2*p.i^2+2*p.i*p.i.ilj+2*p.i.ilj^2)
      il <- li <- ai*al*p.k*p.i.ilj*p.l.lj*((p.j.lj-p.l-p.l.ilj)*(1-p.i-p.i.ilj)+p.i*p.l+p.i.ilj*p.l.ilj)
      ij <- ji <- ai*aj*p.k*p.i.ilj*p.l.lj*((-p.j.lj-p.j-p.j.ilj)*(1-p.i-p.i.ilj)+p.i*p.j+p.i.ilj*p.j.ilj)
      ll <- al^2*p.k*p.i.ilj*p.l.lj*((p.j.lj-p.l-p.l.ilj)^2-p.j.lj*p.l.lj-p.l.ilj*(1-p.l.ilj)-p.l*(1-p.l))
      lj <- jl <- al*aj*p.k*p.i.ilj*p.l.lj*((p.l+p.l.ilj-p.j.lj)*(p.j+p.j.ilj+p.j.lj)+p.l.ilj*p.j.ilj+p.j*p.l+p.l.lj*p.j.lj)
      jj <- aj^2*p.k*p.i.ilj*p.l.lj*((p.j+p.j.ilj+p.j.lj)^2-p.j.ilj*(1-p.j.ilj)-p.j*(1-p.j)-p.j.lj*(1-p.j.lj))
    }else if(identical(x,c(3,2,1,4))){
      kk <- ak^2*p.k*p.j.jil*(1-2*p.k)*p.i.il*(1-p.k)
      kj <- jk <- ak*aj*p.k*p.j.jil*(1-2*p.j)*p.i.il*(1-p.k)
      ki <- ik <- ak*ai*p.k*p.j.jil*p.i.il*(2*p.i*p.k-2*p.i+(1-p.k)*p.l.il)
      kl <- lk <- ak*al*p.k*p.j.jil*p.i.il*(2*p.l*p.k-2*p.l-(1-p.k)*p.l.il)
      jj <- aj^2*p.k*p.j.jil*p.i.il*(1-3*p.j+2*p.j^2+2*p.j*p.j.jil+2*p.j.jil^2)
      ji <- ij <- aj*ai*p.k*p.j.jil*p.i.il*((p.l.il-p.i-p.i.jil)*(1-p.j-p.j.jil)+p.j*p.i+p.j.jil*p.i.jil)
      jl <- lj <- aj*al*p.k*p.j.jil*p.i.il*((-p.l.il-p.l-p.l.jil)*(1-p.j-p.j.jil)+p.j*p.l+p.j.jil*p.l.jil)
      ii <- ai^2*p.k*p.j.jil*p.i.il*((p.l.il-p.i-p.i.jil)^2-p.l.il*p.i.il-p.i.jil*(1-p.i.jil)-p.i*(1-p.i))
      il <- li <- ai*al*p.k*p.j.jil*p.i.il*((p.i+p.i.jil-p.l.il)*(p.l+p.l.jil+p.l.il)+p.i.jil*p.l.jil+p.l*p.i+p.i.il*p.l.il)
      ll <- al^2*p.k*p.j.jil*p.i.il*((p.l+p.l.jil+p.l.il)^2-p.l.jil*(1-p.l.jil)-p.l*(1-p.l)-p.l.il*(1-p.l.il))
    }else if(identical(x,c(3,2,4,1))){
      kk <- ak^2*p.k*p.j.jli*(1-2*p.k)*p.l.li*(1-p.k)
      kj <- jk <- ak*aj*p.k*p.j.jli*(1-2*p.j)*p.l.li*(1-p.k)
      kl <- lk <- ak*al*p.k*p.j.jli*p.l.li*(2*p.l*p.k-2*p.l+(1-p.k)*p.i.li)
      ki <- ik <- ak*ai*p.k*p.j.jli*p.l.li*(2*p.i*p.k-2*p.i-(1-p.k)*p.i.li)
      jj <- aj^2*p.k*p.j.jli*p.l.li*(1-3*p.j+2*p.j^2+2*p.j*p.j.jli+2*p.j.jli^2)
      jl <- lj <- aj*al*p.k*p.j.jli*p.l.li*((p.i.li-p.l-p.l.jli)*(1-p.j-p.j.jli)+p.j*p.l+p.j.jli*p.l.jli)
      ji <- ij <- aj*ai*p.k*p.j.jli*p.l.li*((-p.i.li-p.i-p.i.jli)*(1-p.j-p.j.jli)+p.j*p.i+p.j.jli*p.i.jli)
      ll <- al^2*p.k*p.j.jli*p.l.li*((p.i.li-p.l-p.l.jli)^2-p.i.li*p.l.li-p.l.jli*(1-p.l.jli)-p.l*(1-p.l))
      li <- il <- al*ai*p.k*p.j.jli*p.l.li*((p.l+p.l.jli-p.i.li)*(p.i+p.i.jli+p.i.li)+p.l.jli*p.i.jli+p.i*p.l+p.l.li*p.i.li)
      ii <- ai^2*p.k*p.j.jli*p.l.li*((p.i+p.i.jli+p.i.li)^2-p.i.jli*(1-p.i.jli)-p.i*(1-p.i)-p.i.li*(1-p.i.li))
    }else if(identical(x,c(3,4,1,2))){
      kk <- ak^2*p.k*p.l.lij*(1-2*p.k)*p.i.ij*(1-p.k)
      kl <- lk <- ak*al*p.k*p.l.lij*(1-2*p.l)*p.i.ij*(1-p.k)
      ki <- ik <- ak*ai*p.k*p.l.lij*p.i.ij*(2*p.i*p.k-2*p.i+(1-p.k)*p.j.ij)
      kj <- jk <- ak*aj*p.k*p.l.lij*p.i.ij*(2*p.j*p.k-2*p.j-(1-p.k)*p.j.ij)
      ll <- al^2*p.k*p.l.lij*p.i.ij*(1-3*p.l+2*p.l^2+2*p.l*p.l.lij+2*p.l.lij^2)
      li <- il <- al*ai*p.k*p.l.lij*p.i.ij*((p.j.ij-p.i-p.i.lij)*(1-p.l-p.l.lij)+p.l*p.i+p.l.lij*p.i.lij)
      lj <- jl <- al*aj*p.k*p.l.lij*p.i.ij*((-p.j.ij-p.j-p.j.lij)*(1-p.l-p.l.lij)+p.l*p.j+p.l.lij*p.j.lij)
      ii <- ai^2*p.k*p.l.lij*p.i.ij*((p.j.ij-p.i-p.i.lij)^2-p.j.ij*p.i.ij-p.i.lij*(1-p.i.lij)-p.i*(1-p.i))
      ij <- ji <- ai*aj*p.k*p.l.lij*p.i.ij*((p.i+p.i.lij-p.j.ij)*(p.j+p.j.lij+p.j.ij)+p.i.lij*p.j.lij+p.j*p.i+p.i.ij*p.j.ij)
      jj <- aj^2*p.k*p.l.lij*p.i.ij*((p.j+p.j.lij+p.j.ij)^2-p.j.lij*(1-p.j.lij)-p.j*(1-p.j)-p.j.ij*(1-p.j.ij))
    }else if(identical(x,c(3,4,2,1))){
      kk <- ak^2*p.k*p.l.lji*(1-2*p.k)*p.j.ji*(1-p.k)
      kl <- lk <- ak*al*p.k*p.l.lji*(1-2*p.l)*p.j.ji*(1-p.k)
      kj <- jk <- ak*aj*p.k*p.l.lji*p.j.ji*(2*p.j*p.k-2*p.j+(1-p.k)*p.i.ji)
      ki <- ik <- ak*ai*p.k*p.l.lji*p.j.ji*(2*p.i*p.k-2*p.i-(1-p.k)*p.i.ji)
      ll <- al^2*p.k*p.l.lji*p.j.ji*(1-3*p.l+2*p.l^2+2*p.l*p.l.lji+2*p.l.lji^2)
      lj <- jl <- al*aj*p.k*p.l.lji*p.j.ji*((p.i.ji-p.j-p.j.lji)*(1-p.l-p.l.lji)+p.l*p.j+p.l.lji*p.j.lji)
      li <- il <- al*ai*p.k*p.l.lji*p.j.ji*((-p.i.ji-p.i-p.i.lji)*(1-p.l-p.l.lji)+p.l*p.i+p.l.lji*p.i.lji)
      jj <- aj^2*p.k*p.l.lji*p.j.ji*((p.i.ji-p.j-p.j.lji)^2-p.i.ji*p.j.ji-p.j.lji*(1-p.j.lji)-p.j*(1-p.j))
      ji <- ij <- aj*ai*p.k*p.l.lji*p.j.ji*((p.j+p.j.lji-p.i.ji)*(p.i+p.i.lji+p.i.ji)+p.j.lji*p.i.lji+p.i*p.j+p.j.ji*p.i.ji)
      ii <- ai^2*p.k*p.l.lji*p.j.ji*((p.i+p.i.lji+p.i.ji)^2-p.i.lji*(1-p.i.lji)-p.i*(1-p.i)-p.i.ji*(1-p.i.ji))
    }else if(identical(x,c(4,1,2,3))){
      ll <- al^2*p.l*p.i.ijk*(1-2*p.l)*p.j.jk*(1-p.l)
      li <- il <- al*ai*p.l*p.i.ijk*(1-2*p.i)*p.j.jk*(1-p.l)
      lj <- jl <- al*aj*p.l*p.i.ijk*p.j.jk*(2*p.j*p.l-2*p.j+(1-p.l)*p.k.jk)
      lk <- kl <- al*ak*p.l*p.i.ijk*p.j.jk*(2*p.k*p.l-2*p.k-(1-p.l)*p.k.jk)
      ii <- ai^2*p.l*p.i.ijk*p.j.jk*(1-3*p.i+2*p.i^2+2*p.i*p.i.ijk+2*p.i.ijk^2)
      ij <- ji <- ai*aj*p.l*p.i.ijk*p.j.jk*((p.k.jk-p.j-p.j.ijk)*(1-p.i-p.i.ijk)+p.i*p.j+p.i.ijk*p.j.ijk)
      ik <- ki <- ai*ak*p.l*p.i.ijk*p.j.jk*((-p.k.jk-p.k-p.k.ijk)*(1-p.i-p.i.ijk)+p.i*p.k+p.i.ijk*p.k.ijk)
      jj <- aj^2*p.l*p.i.ijk*p.j.jk*((p.k.jk-p.j-p.j.ijk)^2-p.k.jk*p.j.jk-p.j.ijk*(1-p.j.ijk)-p.j*(1-p.j))
      jk <- kj <- aj*ak*p.l*p.i.ijk*p.j.jk*((p.j+p.j.ijk-p.k.jk)*(p.k+p.k.ijk+p.k.jk)+p.j.ijk*p.k.ijk+p.k*p.j+p.j.jk*p.k.jk)
      kk <- ak^2*p.l*p.i.ijk*p.j.jk*((p.k+p.k.ijk+p.k.jk)^2-p.k.ijk*(1-p.k.ijk)-p.k*(1-p.k)-p.k.jk*(1-p.k.jk))
    }else if(identical(x,c(4,1,3,2))){
      ll <- al^2*p.l*p.i.ikj*(1-2*p.l)*p.k.kj*(1-p.l)
      li <- il <- al*ai*p.l*p.i.ikj*(1-2*p.i)*p.k.kj*(1-p.l)
      lk <- kl <- al*ak*p.l*p.i.ikj*p.k.kj*(2*p.k*p.l-2*p.k+(1-p.l)*p.j.kj)
      lj <- jl <- al*aj*p.l*p.i.ikj*p.k.kj*(2*p.j*p.l-2*p.j-(1-p.l)*p.j.kj)
      ii <- ai^2*p.l*p.i.ikj*p.k.kj*(1-3*p.i+2*p.i^2+2*p.i*p.i.ikj+2*p.i.ikj^2)
      ik <- ki <- ai*ak*p.l*p.i.ikj*p.k.kj*((p.j.kj-p.k-p.k.ikj)*(1-p.i-p.i.ikj)+p.i*p.k+p.i.ikj*p.k.ikj)
      ij <- ji <- ai*aj*p.l*p.i.ikj*p.k.kj*((-p.j.kj-p.j-p.j.ikj)*(1-p.i-p.i.ikj)+p.i*p.j+p.i.ikj*p.j.ikj)
      kk <- ak^2*p.l*p.i.ikj*p.k.kj*((p.j.kj-p.k-p.k.ikj)^2-p.j.kj*p.k.kj-p.k.ikj*(1-p.k.ikj)-p.k*(1-p.k))
      kj <- jk <- ak*aj*p.l*p.i.ikj*p.k.kj*((p.k+p.k.ikj-p.j.kj)*(p.j+p.j.ikj+p.j.kj)+p.k.ikj*p.j.ikj+p.j*p.k+p.k.kj*p.j.kj)
      jj <- aj^2*p.l*p.i.ikj*p.k.kj*((p.j+p.j.ikj+p.j.kj)^2-p.j.ikj*(1-p.j.ikj)-p.j*(1-p.j)-p.j.kj*(1-p.j.kj))
    }else if(identical(x,c(4,2,1,3))){
      ll <- al^2*p.l*p.j.jik*(1-2*p.l)*p.i.ik*(1-p.l)
      lj <- jl <- al*aj*p.l*p.j.jik*(1-2*p.j)*p.i.ik*(1-p.l)
      li <- il <- al*ai*p.l*p.j.jik*p.i.ik*(2*p.i*p.l-2*p.i+(1-p.l)*p.k.ik)
      lk <- kl <- al*ak*p.l*p.j.jik*p.i.ik*(2*p.k*p.l-2*p.k-(1-p.l)*p.k.ik)
      jj <- aj^2*p.l*p.j.jik*p.i.ik*(1-3*p.j+2*p.j^2+2*p.j*p.j.jik+2*p.j.jik^2)
      ji <- ij <- aj*ai*p.l*p.j.jik*p.i.ik*((p.k.ik-p.i-p.i.jik)*(1-p.j-p.j.jik)+p.j*p.i+p.j.jik*p.i.jik)
      jk <- kj <- aj*ak*p.l*p.j.jik*p.i.ik*((-p.k.ik-p.k-p.k.jik)*(1-p.j-p.j.jik)+p.j*p.k+p.j.jik*p.k.jik)
      ii <- ai^2*p.l*p.j.jik*p.i.ik*((p.k.ik-p.i-p.i.jik)^2-p.k.ik*p.i.ik-p.i.jik*(1-p.i.jik)-p.i*(1-p.i))
      ik <- ki <- ai*ak*p.l*p.j.jik*p.i.ik*((p.i+p.i.jik-p.k.ik)*(p.k+p.k.jik+p.k.ik)+p.i.jik*p.k.jik+p.k*p.i+p.i.ik*p.k.ik)
      kk <- ak^2*p.l*p.j.jik*p.i.ik*((p.k+p.k.jik+p.k.ik)^2-p.k.jik*(1-p.k.jik)-p.k*(1-p.k)-p.k.ik*(1-p.k.ik))
    }else if(identical(x,c(4,2,3,1))){
      ll <- al^2*p.l*p.j.jki*(1-2*p.l)*p.k.ki*(1-p.l)
      lj <- jl <- al*aj*p.l*p.j.jki*(1-2*p.j)*p.k.ki*(1-p.l)
      lk <- kl <- al*ak*p.l*p.j.jki*p.k.ki*(2*p.k*p.l-2*p.k+(1-p.l)*p.i.ki)
      li <- il <- al*ai*p.l*p.j.jki*p.k.ki*(2*p.i*p.l-2*p.i-(1-p.l)*p.i.ki)
      jj <- aj^2*p.l*p.j.jki*p.k.ki*(1-3*p.j+2*p.j^2+2*p.j*p.j.jki+2*p.j.jki^2)
      jk <- kj <- aj*ak*p.l*p.j.jki*p.k.ki*((p.i.ki-p.k-p.k.jki)*(1-p.j-p.j.jki)+p.j*p.k+p.j.jki*p.k.jki)
      ji <- ij <- aj*ai*p.l*p.j.jki*p.k.ki*((-p.i.ki-p.i-p.i.jki)*(1-p.j-p.j.jki)+p.j*p.i+p.j.jki*p.i.jki)
      kk <- ak^2*p.l*p.j.jki*p.k.ki*((p.i.ki-p.k-p.k.jki)^2-p.i.ki*p.k.ki-p.k.jki*(1-p.k.jki)-p.k*(1-p.k))
      ki <- ik <- ak*ai*p.l*p.j.jki*p.k.ki*((p.k+p.k.jki-p.i.ki)*(p.i+p.i.jki+p.i.ki)+p.k.jki*p.i.jki+p.i*p.k+p.k.ki*p.i.ki)
      ii <- ai^2*p.l*p.j.jki*p.k.ki*((p.i+p.i.jki+p.i.ki)^2-p.i.jki*(1-p.i.jki)-p.i*(1-p.i)-p.i.ki*(1-p.i.ki))
    }else if(identical(x,c(4,3,1,2))){
      ll <- al^2*p.l*p.k.kij*(1-2*p.l)*p.i.ij*(1-p.l)
      lk <- kl <- al*ak*p.l*p.k.kij*(1-2*p.k)*p.i.ij*(1-p.l)
      li <- il <- al*ai*p.l*p.k.kij*p.i.ij*(2*p.i*p.l-2*p.i+(1-p.l)*p.j.ij)
      lj <- jl <- al*aj*p.l*p.k.kij*p.i.ij*(2*p.j*p.l-2*p.j-(1-p.l)*p.j.ij)
      kk <- ak^2*p.l*p.k.kij*p.i.ij*(1-3*p.k+2*p.k^2+2*p.k*p.k.kij+2*p.k.kij^2)
      ki <- ik <- ak*ai*p.l*p.k.kij*p.i.ij*((p.j.ij-p.i-p.i.kij)*(1-p.k-p.k.kij)+p.k*p.i+p.k.kij*p.i.kij)
      kj <- jk <- ak*aj*p.l*p.k.kij*p.i.ij*((-p.j.ij-p.j-p.j.kij)*(1-p.k-p.k.kij)+p.k*p.j+p.k.kij*p.j.kij)
      ii <- ai^2*p.l*p.k.kij*p.i.ij*((p.j.ij-p.i-p.i.kij)^2-p.j.ij*p.i.ij-p.i.kij*(1-p.i.kij)-p.i*(1-p.i))
      ij <- ji <- ai*aj*p.l*p.k.kij*p.i.ij*((p.i+p.i.kij-p.j.ij)*(p.j+p.j.kij+p.j.ij)+p.i.kij*p.j.kij+p.j*p.i+p.i.ij*p.j.ij)
      jj <- aj^2*p.l*p.k.kij*p.i.ij*((p.j+p.j.kij+p.j.ij)^2-p.j.kij*(1-p.j.kij)-p.j*(1-p.j)-p.j.ij*(1-p.j.ij))
    }else if(identical(x,c(4,3,2,1))){
      ll <- al^2*p.l*p.k.kji*(1-2*p.l)*p.j.ji*(1-p.l)
      lk <- kl <- al*ak*p.l*p.k.kji*(1-2*p.k)*p.j.ji*(1-p.l)
      lj <- jl <- al*aj*p.l*p.k.kji*p.j.ji*(2*p.j*p.l-2*p.j+(1-p.l)*p.i.ji)
      li <- il <- al*ai*p.l*p.k.kji*p.j.ji*(2*p.i*p.l-2*p.i-(1-p.l)*p.i.ji)
      kk <- ak^2*p.l*p.k.kji*p.j.ji*(1-3*p.k+2*p.k^2+2*p.k*p.k.kji+2*p.k.kji^2)
      kj <- jk <- ak*aj*p.l*p.k.kji*p.j.ji*((p.i.ji-p.j-p.j.kji)*(1-p.k-p.k.kji)+p.k*p.j+p.k.kji*p.j.kji)
      ki <- ik <- ak*ai*p.l*p.k.kji*p.j.ji*((-p.i.ji-p.i-p.i.kji)*(1-p.k-p.k.kji)+p.k*p.i+p.k.kji*p.i.kji)
      jj <- aj^2*p.l*p.k.kji*p.j.ji*((p.i.ji-p.j-p.j.kji)^2-p.i.ji*p.j.ji-p.j.kji*(1-p.j.kji)-p.j*(1-p.j))
      ji <- ij <- aj*ai*p.l*p.k.kji*p.j.ji*((p.j+p.j.kji-p.i.ji)*(p.i+p.i.kji+p.i.ji)+p.j.kji*p.i.kji+p.i*p.j+p.j.ji*p.i.ji)
      ii <- ai^2*p.l*p.k.kji*p.j.ji*((p.i+p.i.kji+p.i.ji)^2-p.i.kji*(1-p.i.kji)-p.i*(1-p.i)-p.i.ji*(1-p.i.ji))
    }
  h <- diag(c(ii,jj,kk))
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


Ib_theta_2PL_4_rank <- function(theta,a,d){
  all.x <- list(c(1,2,3,4),c(1,2,4,3),c(1,3,2,4),c(1,3,4,2),c(1,4,2,3),c(1,4,3,2),
                c(2,1,3,4),c(2,1,4,3),c(2,3,1,4),c(2,3,4,1),c(2,4,1,3),c(2,4,3,1),
                c(3,1,2,4),c(3,1,4,2),c(3,2,1,4),c(3,2,4,1),c(3,4,1,2),c(3,4,2,1),
                c(4,1,2,3),c(4,1,3,2),c(4,2,1,3),c(4,2,3,1),c(4,3,1,2),c(4,3,2,1))
  ret <- 0
  for(x in 1:length(all.x)){
    g <- grad_theta_2PL_4_rank(theta,a,d,x=all.x[[x]])
    ret <- ret + tcrossprod(g,g)/px_2PL_4_rank(theta,a,d,xx=all.x[[x]])-hes_theta_2PL_4_rank(theta,a,d,x=all.x[[x]])
  }
  ret
}


test.info_theta_2PL_4_rank <- function(theta,item.par,BID){
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
    complete.info.b[dim.b,dim.b] <- info[[b]] <- Ib_theta_2PL_4_rank(theta.b,a.b,d.b)
    I <- I+complete.info.b
    comp.info[[b]] <- complete.info.b
  }
  se <- sqrt(diag(solve(Reduce("+",comp.info))))
  list(comp.info=comp.info,se =se, info.b=info,I=I)
}

