px_2PL_4_mole <- function(theta,a,d,x=c(1,2,3,4)){ #a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)
  z <- exp(a*theta-d)
  prod(z[x[-4]])/(sum(z)*sum(z[x[-1]])*sum(z[x[-1:-2]]))
}

grad_theta_2PL_4_mole <- function(theta,a,d,x=c(1,2,3,4)){
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
    g.i <- ai*p.i*p.j*p.k.kl
    g.j <- aj*p.i*p.j.jkl*(p.k.jkl-p.j*p.k.kl)
    g.k <- ak*p.i*p.j.jkl*p.k.kl*(p.k.kl-p.k.jkl-p.k)
    g.l <- -al*p.i*p.j.jkl*p.k.kl*(p.l+p.l.jkl+p.l.kl)
  }else if(identical(x,c(1,2,4,3))){
    g.i <- ai*p.i*p.j*p.l.kl
    g.j <- aj*p.i*p.j.jkl*(p.l.jkl-p.j*p.l.kl)
    g.l <- al*p.i*p.j.jkl*p.l.kl*(p.l.kl-p.l.jkl-p.l)
    g.k <- -ak*p.i*p.j.jkl*p.l.kl*(p.k+p.k.jkl+p.k.kl)
  }else if(identical(x,c(1,3,2,4))){
    g.i <- ai*p.i*p.k*p.j.jl
    g.k <- ak*p.i*p.k.jkl*(p.j.jkl-p.k*p.j.jl)
    g.j <- aj*p.i*p.k.jkl*p.j.jl*(p.j.jl-p.j.jkl-p.j)
    g.l <- -al*p.i*p.k.jkl*p.j.jl*(p.l+p.l.jkl+p.l.jl)
  }else if(identical(x,c(1,3,4,2))){
    g.i <- ai*p.i*p.k*p.l.jl
    g.k <- ak*p.i*p.k.klj*(p.l.klj-p.k*p.l.lj)
    g.l <- al*p.i*p.k.klj*p.l.lj*(p.l.lj-p.l.klj-p.l)
    g.j <- -aj*p.i*p.k.klj*p.l.lj*(p.j+p.j.klj+p.j.lj)
  }else if(identical(x,c(1,4,2,3))){
    g.i <- ai*p.i*p.l*p.j.jk
    g.l <- al*p.i*p.l.ljk*(p.j.ljk-p.l*p.j.jk)
    g.j <- aj*p.i*p.l.ljk*p.j.jk*(p.j.jk-p.j.ljk-p.j)
    g.k <- -ak*p.i*p.l.ljk*p.j.jk*(p.k+p.k.ljk+p.k.jk)
  }else if(identical(x,c(1,4,3,2))){
    g.i <- ai*p.i*p.l*p.k.kj
    g.l <- al*p.i*p.l.lkj*(p.k.lkj-p.l*p.k.kj)
    g.k <- ak*p.i*p.l.lkj*p.k.kj*(p.k.kj-p.k.lkj-p.k)
    g.j <- -aj*p.i*p.l.lkj*p.k.kj*(p.j+p.j.lkj+p.j.kj)
  }else if(identical(x,c(2,1,3,4))){
    g.j <- aj*p.j*p.i*p.k.kl
    g.i <- ai*p.j*p.i.ikl*(p.k.ikl-p.i*p.k.kl)
    g.k <- ak*p.j*p.i.ikl*p.k.kl*(p.k.kl-p.k.ikl-p.k)
    g.l <- -al*p.j*p.i.ikl*p.k.kl*(p.l+p.l.ikl+p.l.kl)
  }else if(identical(x,c(2,1,4,3))){
    g.j <- aj*p.j*p.i*p.l.lk
    g.i <- ai*p.j*p.i.ilk*(p.l.ilk-p.i*p.l.lk)
    g.l <- al*p.j*p.i.ilk*p.l.lk*(p.l.lk-p.l.ilk-p.l)
    g.k <- -ak*p.j*p.i.ilk*p.l.lk*(p.k+p.k.ilk+p.k.lk)
  }else if(identical(x,c(2,3,1,4))){
    g.j <- aj*p.j*p.k*p.i.il
    g.k <- ak*p.j*p.k.kil*(p.i.kil-p.k*p.i.il)
    g.i <- ai*p.j*p.k.kil*p.i.il*(p.i.il-p.i.kil-p.i)
    g.l <- -al*p.j*p.k.kil*p.i.il*(p.l+p.l.kil+p.l.il)
  }else if(identical(x,c(2,3,4,1))){
    g.j <- aj*p.j*p.k*p.l.li
    g.k <- ak*p.j*p.k.kli*(p.l.kli-p.k*p.l.li)
    g.l <- al*p.j*p.k.kli*p.l.li*(p.l.li-p.l.kli-p.l)
    g.i <- -ai*p.j*p.k.kli*p.l.li*(p.i+p.i.kli+p.i.li)
  }else if(identical(x,c(2,4,1,3))){
    g.j <- aj*p.j*p.l*p.i.ik
    g.l <- al*p.j*p.l.lik*(p.i.lik-p.l*p.i.ik)
    g.i <- ai*p.j*p.l.lik*p.i.ik*(p.i.ik-p.i.lik-p.i)
    g.k <- -ak*p.j*p.l.lik*p.i.ik*(p.k+p.k.lik+p.k.ik)
  }else if(identical(x,c(2,4,3,1))){
    g.i <- ai*p.i*p.j*p.k.kl
    g.j <- aj*p.i*p.j.jkl*(p.k.jkl-p.j*p.k.kl)
    g.k <- ak*p.i*p.j.jkl*p.k.kl*(p.k.kl-p.k.jkl-p.k)
    g.l <- -al*p.i*p.j.jkl*p.k.kl*(p.l+p.l.jkl+p.l.kl)
  }else if(identical(x,c(3,1,2,4))){
    g.k <- ak*p.k*p.i*p.j.jl
    g.i <- ai*p.k*p.i.ijl*(p.j.ijl-p.i*p.j.jl)
    g.j <- aj*p.k*p.i.ijl*p.j.jl*(p.j.jl-p.j.ijl-p.j)
    g.l <- -al*p.k*p.i.ijl*p.j.jl*(p.l+p.l.ijl+p.l.jl)
  }else if(identical(x,c(3,1,4,2))){
    g.k <- ak*p.k*p.i*p.l.lj
    g.i <- ai*p.k*p.i.ilj*(p.l.ilj-p.i*p.l.lj)
    g.l <- al*p.k*p.i.ilj*p.l.lj*(p.l.lj-p.l.ilj-p.l)
    g.j <- -aj*p.k*p.i.ilj*p.l.lj*(p.j+p.j.ilj+p.j.lj)
  }else if(identical(x,c(3,2,1,4))){
    g.k <- ak*p.k*p.j*p.i.il
    g.j <- aj*p.k*p.j.jil*(p.i.jil-p.j*p.i.il)
    g.i <- ai*p.k*p.j.jil*p.i.il*(p.i.il-p.i.jil-p.i)
    g.l <- -al*p.k*p.j.jil*p.i.il*(p.l+p.l.jil+p.l.il)
  }else if(identical(x,c(3,2,4,1))){
    g.k <- ak*p.k*p.j*p.l.li
    g.j <- aj*p.k*p.j.jli*(p.l.jli-p.j*p.l.li)
    g.l <- al*p.k*p.j.jli*p.l.li*(p.l.li-p.l.jli-p.l)
    g.i <- -ai*p.k*p.j.jli*p.l.li*(p.i+p.i.jli+p.i.li)
  }else if(identical(x,c(3,4,1,2))){
    g.k <- ak*p.k*p.l*p.i.ij
    g.l <- al*p.k*p.l.lij*(p.i.lij-p.l*p.i.ij)
    g.i <- ai*p.k*p.l.lij*p.i.ij*(p.i.ij-p.i.lij-p.i)
    g.j <- -aj*p.k*p.l.lij*p.i.ij*(p.j+p.j.lij+p.j.ij)
  }else if(identical(x,c(3,4,2,1))){
    g.k <- ak*p.k*p.l*p.j.ji
    g.l <- al*p.k*p.l.lji*(p.j.lji-p.l*p.j.ji)
    g.j <- aj*p.k*p.l.lji*p.j.ji*(p.j.ji-p.j.lji-p.j)
    g.i <- -ai*p.k*p.l.lji*p.j.ji*(p.i+p.i.lji+p.i.ji)
  }else if(identical(x,c(4,1,2,3))){
    g.l <- al*p.l*p.i*p.j.jk
    g.i <- ai*p.l*p.i.ijk*(p.j.ijk-p.i*p.j.jk)
    g.j <- aj*p.l*p.i.ijk*p.j.jk*(p.j.jk-p.j.ijk-p.j)
    g.k <- -ak*p.l*p.i.ijk*p.j.jk*(p.k+p.k.ijk+p.k.jk)
  }else if(identical(x,c(4,1,3,2))){
    g.l <- al*p.l*p.i*p.k.kj
    g.i <- ai*p.l*p.i.ikj*(p.k.ikj-p.i*p.k.kj)
    g.k <- ak*p.l*p.i.ikj*p.k.kj*(p.k.kj-p.k.ikj-p.k)
    g.j <- -aj*p.l*p.i.ikj*p.k.kj*(p.j+p.j.ikj+p.j.kj)
  }else if(identical(x,c(4,2,1,3))){
    g.l <- al*p.l*p.j*p.i.ik
    g.j <- aj*p.l*p.j.jik*(p.i.jik-p.j*p.i.ik)
    g.i <- ai*p.l*p.j.jik*p.i.ik*(p.i.ik-p.i.jik-p.i)
    g.k <- -ak*p.l*p.j.jik*p.i.ik*(p.k+p.k.jik+p.k.ik)
  }else if(identical(x,c(4,2,3,1))){
    g.l <- al*p.l*p.j*p.k.ki
    g.j <- aj*p.l*p.j.jki*(p.k.jki-p.j*p.k.ki)
    g.k <- ak*p.l*p.j.jki*p.k.ki*(p.k.ki-p.k.jki-p.k)
    g.i <- -ai*p.l*p.j.jki*p.k.ki*(p.i+p.i.jki+p.i.ki)
  }else if(identical(x,c(4,3,1,2))){
    g.l <- al*p.l*p.k*p.i.ij
    g.k <- ak*p.l*p.k.kij*(p.i.kij-p.k*p.i.ij)
    g.i <- ai*p.l*p.k.kij*p.i.ij*(p.i.ij-p.i.kij-p.i)
    g.j <- -aj*p.l*p.k.kij*p.i.ij*(p.j+p.j.kij+p.j.ij)
  }else if(identical(x,c(4,3,2,1))){
    g.l <- al*p.l*p.k*p.j.ji
    g.k <- ak*p.l*p.k.kji*(p.j.kji-p.k*p.j.ji)
    g.j <- aj*p.l*p.k.kji*p.j.ji*(p.j.ji-p.j.kji-p.j)
    g.i <- -ai*p.l*p.k.kji*p.j.ji*(p.i+p.i.kji+p.i.ji)
  }

  c(g.i,g.j,g.k,g.l)
}

hes_theta_2PL_4_mole <- function(theta,a,d,x=c(1,2,3,4)){
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
    ii <- ai^2*p.i*p.j*(1-2*p.i)*p.k.kl
    ij <- ji <- ai*aj*p.i*p.j*(1-2*p.j)*p.k.kl
    ik <- ki <- ai*ak*p.i*p.j*p.k.kl*(p.k.kl-2*p.k)
    il <- li <- -ai*al*p.i*p.j*p.k.kl*(p.l.kl+2*p.l)
    jj <- aj^2*p.i*p.j.jkl*((2*p.j-2+p.j.jkl)*p.j*p.k.kl+(1-p.j)*pk.jkl-2*p.k.jkl*p.j.jkl)
    jk <- kj <- aj*ak*p.i*p.j.jkl*(p.k.jkl*(1-p.k)+2*p.k*p.j*p.k.kl-2*p.k.jkl^2+(p.k.jkl-p.l.kl)*p.j*p.k.kl)
    jl <- lj <- aj*al*p.i*p.j.jkl*(2*p.l*p.j*p.k.kl-p.l*p.k.jkl-2*p.l.jkl*p.k.jkl+p.j*p.k.kl*(p.l.jkl+p.l.kl))
    kk <- 2*ak^2*p.i*p.j.jkl*p.k.kl*((p.k-1)*(p.k+p.k.jkl)+p.k.jkl^2+p.l.kl*p.k.kl)
    kl <- lk <- ak*al*p.i*p.j.jkl*p.k.kl*((p.k+p.k.jkl-p.k.kl)*(p.l+p.l.jkl+p.l.kl)+p.k.jkl*p.l.jkl+p.l*p.k-p.k.kl*p.l.kl)
    ll <- al^2*p.i*p.j.jkl*p.k.kl*((p.l+p.l.jkl+p.l.kl)^2-p.l.jkl*(1-p.l.jkl)-p.l*(1-p.l)-p.l.kl*(1-p.l.kl))
  }else if(identical(x,c(1,2,4,3))){
    ii <- ai^2*p.i*p.j*(1-2*p.i)*p.l.lk
    ij <- ji <- ai*aj*p.i*p.j*(1-2*p.j)*p.l.lk
    il <- li <- ai*al*p.i*p.j*p.l.lk*(p.l.lk-2*p.l)
    ik <- ki <- -ai*ak*p.i*p.j*p.l.lk*(p.k.lk+2*p.k)
    jj <- aj^2*p.i*p.j.jlk*((2*p.j-2+p.j.jlk)*p.j*p.l.lk+(1-p.j)*pl.jlk-2*p.l.jlk*p.j.jlk)
    jl <- lj <- aj*al*p.i*p.j.jlk*(p.l.jlk*(1-p.l)+2*p.l*p.j*p.l.lk-2*p.l.jlk^2+(p.l.jlk-p.k.lk)*p.j*p.l.lk)
    jk <- kj <- aj*ak*p.i*p.j.jlk*(2*p.k*p.j*p.l.lk-p.k*p.l.jlk-2*p.k.jlk*p.l.jlk+p.j*p.l.lk*(p.k.jlk+p.k.lk))
    ll <- 2*al^2*p.i*p.j.jlk*p.l.lk*((p.l-1)*(p.l+p.l.jlk)+p.l.jlk^2+p.k.lk*p.l.lk)
    lk <- kl <- al*ak*p.i*p.j.jlk*p.l.lk*((p.l+p.l.jlk-p.l.lk)*(p.k+p.k.jlk+p.k.lk)+p.l.jlk*p.k.jlk+p.k*p.l-p.l.lk*p.k.lk)
    kk <- ak^2*p.i*p.j.jlk*p.l.lk*((p.k+p.k.jlk+p.k.lk)^2-p.k.jlk*(1-p.k.jlk)-p.k*(1-p.k)-p.k.lk*(1-p.k.lk))
  }else if(identical(x,c(1,3,2,4))){
    ii <- ai^2*p.i*p.k*(1-2*p.i)*p.j.jl
    ik <- ki <- ai*ak*p.i*p.k*(1-2*p.k)*p.j.jl
    ij <- ji <- ai*aj*p.i*p.k*p.j.jl*(p.j.jl-2*p.j)
    il <- li <- -ai*al*p.i*p.k*p.j.jl*(p.l.jl+2*p.l)
    kk <- ak^2*p.i*p.k.kjl*((2*p.k-2+p.k.kjl)*p.k*p.j.jl+(1-p.k)*pj.kjl-2*p.j.kjl*p.k.kjl)
    kj <- jk <- ak*aj*p.i*p.k.kjl*(p.j.kjl*(1-p.j)+2*p.j*p.k*p.j.jl-2*p.j.kjl^2+(p.j.kjl-p.l.jl)*p.k*p.j.jl)
    kl <- lk <- ak*al*p.i*p.k.kjl*(2*p.l*p.k*p.j.jl-p.l*p.j.kjl-2*p.l.kjl*p.j.kjl+p.k*p.j.jl*(p.l.kjl+p.l.jl))
    jj <- 2*aj^2*p.i*p.k.kjl*p.j.jl*((p.j-1)*(p.j+p.j.kjl)+p.j.kjl^2+p.l.jl*p.j.jl)
    jl <- lj <- aj*al*p.i*p.k.kjl*p.j.jl*((p.j+p.j.kjl-p.j.jl)*(p.l+p.l.kjl+p.l.jl)+p.j.kjl*p.l.kjl+p.l*p.j-p.j.jl*p.l.jl)
    ll <- al^2*p.i*p.k.kjl*p.j.jl*((p.l+p.l.kjl+p.l.jl)^2-p.l.kjl*(1-p.l.kjl)-p.l*(1-p.l)-p.l.jl*(1-p.l.jl))
  }else if(identical(x,c(1,3,4,2))){
    ii <- ai^2*p.i*p.k*(1-2*p.i)*p.l.lj
    ik <- ki <- ai*ak*p.i*p.k*(1-2*p.k)*p.l.lj
    il <- li <- ai*al*p.i*p.k*p.l.lj*(p.l.lj-2*p.l)
    ij <- ji <- -ai*aj*p.i*p.k*p.l.lj*(p.j.lj+2*p.j)
    kk <- ak^2*p.i*p.k.klj*((2*p.k-2+p.k.klj)*p.k*p.l.lj+(1-p.k)*pl.klj-2*p.l.klj*p.k.klj)
    kl <- lk <- ak*al*p.i*p.k.klj*(p.l.klj*(1-p.l)+2*p.l*p.k*p.l.lj-2*p.l.klj^2+(p.l.klj-p.j.lj)*p.k*p.l.lj)
    kj <- jk <- ak*aj*p.i*p.k.klj*(2*p.j*p.k*p.l.lj-p.j*p.l.klj-2*p.j.klj*p.l.klj+p.k*p.l.lj*(p.j.klj+p.j.lj))
    ll <- 2*al^2*p.i*p.k.klj*p.l.lj*((p.l-1)*(p.l+p.l.klj)+p.l.klj^2+p.j.lj*p.l.lj)
    lj <- jl <- al*aj*p.i*p.k.klj*p.l.lj*((p.l+p.l.klj-p.l.lj)*(p.j+p.j.klj+p.j.lj)+p.l.klj*p.j.klj+p.j*p.l-p.l.lj*p.j.lj)
    jj <- aj^2*p.i*p.k.klj*p.l.lj*((p.j+p.j.klj+p.j.lj)^2-p.j.klj*(1-p.j.klj)-p.j*(1-p.j)-p.j.lj*(1-p.j.lj))
  }else if(identical(x,c(1,4,2,3))){
    ii <- ai^2*p.i*p.l*(1-2*p.i)*p.j.jk
    il <- li <- ai*al*p.i*p.l*(1-2*p.l)*p.j.jk
    ij <- ji <- ai*aj*p.i*p.l*p.j.jk*(p.j.jk-2*p.j)
    ik <- ki <- -ai*ak*p.i*p.l*p.j.jk*(p.k.jk+2*p.k)
    ll <- al^2*p.i*p.l.ljk*((2*p.l-2+p.l.ljk)*p.l*p.j.jk+(1-p.l)*pj.ljk-2*p.j.ljk*p.l.ljk)
    lj <- jl <- al*aj*p.i*p.l.ljk*(p.j.ljk*(1-p.j)+2*p.j*p.l*p.j.jk-2*p.j.ljk^2+(p.j.ljk-p.k.jk)*p.l*p.j.jk)
    lk <- kl <- al*ak*p.i*p.l.ljk*(2*p.k*p.l*p.j.jk-p.k*p.j.ljk-2*p.k.ljk*p.j.ljk+p.l*p.j.jk*(p.k.ljk+p.k.jk))
    jj <- 2*aj^2*p.i*p.l.ljk*p.j.jk*((p.j-1)*(p.j+p.j.ljk)+p.j.ljk^2+p.k.jk*p.j.jk)
    jk <- kj <- aj*ak*p.i*p.l.ljk*p.j.jk*((p.j+p.j.ljk-p.j.jk)*(p.k+p.k.ljk+p.k.jk)+p.j.ljk*p.k.ljk+p.k*p.j-p.j.jk*p.k.jk)
    kk <- ak^2*p.i*p.l.ljk*p.j.jk*((p.k+p.k.ljk+p.k.jk)^2-p.k.ljk*(1-p.k.ljk)-p.k*(1-p.k)-p.k.jk*(1-p.k.jk))
  }else if(identical(x,c(1,4,3,2))){
    ii <- ai^2*p.i*p.l*(1-2*p.i)*p.k.kj
    il <- li <- ai*al*p.i*p.l*(1-2*p.l)*p.k.kj
    ik <- ki <- ai*ak*p.i*p.l*p.k.kj*(p.k.kj-2*p.k)
    ij <- ji <- -ai*aj*p.i*p.l*p.k.kj*(p.j.kj+2*p.j)
    ll <- al^2*p.i*p.l.lkj*((2*p.l-2+p.l.lkj)*p.l*p.k.kj+(1-p.l)*pk.lkj-2*p.k.lkj*p.l.lkj)
    lk <- kl <- al*ak*p.i*p.l.lkj*(p.k.lkj*(1-p.k)+2*p.k*p.l*p.k.kj-2*p.k.lkj^2+(p.k.lkj-p.j.kj)*p.l*p.k.kj)
    lj <- jl <- al*aj*p.i*p.l.lkj*(2*p.j*p.l*p.k.kj-p.j*p.k.lkj-2*p.j.lkj*p.k.lkj+p.l*p.k.kj*(p.j.lkj+p.j.kj))
    kk <- 2*ak^2*p.i*p.l.lkj*p.k.kj*((p.k-1)*(p.k+p.k.lkj)+p.k.lkj^2+p.j.kj*p.k.kj)
    kj <- jk <- ak*aj*p.i*p.l.lkj*p.k.kj*((p.k+p.k.lkj-p.k.kj)*(p.j+p.j.lkj+p.j.kj)+p.k.lkj*p.j.lkj+p.j*p.k-p.k.kj*p.j.kj)
    jj <- aj^2*p.i*p.l.lkj*p.k.kj*((p.j+p.j.lkj+p.j.kj)^2-p.j.lkj*(1-p.j.lkj)-p.j*(1-p.j)-p.j.kj*(1-p.j.kj))
  }else if(identical(x,c(2,1,3,4))){
    jj <- aj^2*p.j*p.i*(1-2*p.j)*p.k.kl
    ji <- ij <- aj*ai*p.j*p.i*(1-2*p.i)*p.k.kl
    jk <- kj <- aj*ak*p.j*p.i*p.k.kl*(p.k.kl-2*p.k)
    jl <- lj <- -aj*al*p.j*p.i*p.k.kl*(p.l.kl+2*p.l)
    ii <- ai^2*p.j*p.i.ikl*((2*p.i-2+p.i.ikl)*p.i*p.k.kl+(1-p.i)*pk.ikl-2*p.k.ikl*p.i.ikl)
    ik <- ki <- ai*ak*p.j*p.i.ikl*(p.k.ikl*(1-p.k)+2*p.k*p.i*p.k.kl-2*p.k.ikl^2+(p.k.ikl-p.l.kl)*p.i*p.k.kl)
    il <- li <- ai*al*p.j*p.i.ikl*(2*p.l*p.i*p.k.kl-p.l*p.k.ikl-2*p.l.ikl*p.k.ikl+p.i*p.k.kl*(p.l.ikl+p.l.kl))
    kk <- 2*ak^2*p.j*p.i.ikl*p.k.kl*((p.k-1)*(p.k+p.k.ikl)+p.k.ikl^2+p.l.kl*p.k.kl)
    kl <- lk <- ak*al*p.j*p.i.ikl*p.k.kl*((p.k+p.k.ikl-p.k.kl)*(p.l+p.l.ikl+p.l.kl)+p.k.ikl*p.l.ikl+p.l*p.k-p.k.kl*p.l.kl)
    ll <- al^2*p.j*p.i.ikl*p.k.kl*((p.l+p.l.ikl+p.l.kl)^2-p.l.ikl*(1-p.l.ikl)-p.l*(1-p.l)-p.l.kl*(1-p.l.kl))
  }else if(identical(x,c(2,1,4,3))){
    jj <- aj^2*p.j*p.i*(1-2*p.j)*p.l.lk
    ji <- ij <- aj*ai*p.j*p.i*(1-2*p.i)*p.l.lk
    jl <- lj <- aj*al*p.j*p.i*p.l.lk*(p.l.lk-2*p.l)
    jk <- kj <- -aj*ak*p.j*p.i*p.l.lk*(p.k.lk+2*p.k)
    ii <- ai^2*p.j*p.i.ilk*((2*p.i-2+p.i.ilk)*p.i*p.l.lk+(1-p.i)*pl.ilk-2*p.l.ilk*p.i.ilk)
    il <- li <- ai*al*p.j*p.i.ilk*(p.l.ilk*(1-p.l)+2*p.l*p.i*p.l.lk-2*p.l.ilk^2+(p.l.ilk-p.k.lk)*p.i*p.l.lk)
    ik <- ki <- ai*ak*p.j*p.i.ilk*(2*p.k*p.i*p.l.lk-p.k*p.l.ilk-2*p.k.ilk*p.l.ilk+p.i*p.l.lk*(p.k.ilk+p.k.lk))
    ll <- 2*al^2*p.j*p.i.ilk*p.l.lk*((p.l-1)*(p.l+p.l.ilk)+p.l.ilk^2+p.k.lk*p.l.lk)
    lk <- kl <- al*ak*p.j*p.i.ilk*p.l.lk*((p.l+p.l.ilk-p.l.lk)*(p.k+p.k.ilk+p.k.lk)+p.l.ilk*p.k.ilk+p.k*p.l-p.l.lk*p.k.lk)
    kk <- ak^2*p.j*p.i.ilk*p.l.lk*((p.k+p.k.ilk+p.k.lk)^2-p.k.ilk*(1-p.k.ilk)-p.k*(1-p.k)-p.k.lk*(1-p.k.lk))
  }else if(identical(x,c(2,3,1,4))){
    jj <- aj^2*p.j*p.k*(1-2*p.j)*p.i.il
    jk <- kj <- aj*ak*p.j*p.k*(1-2*p.k)*p.i.il
    ji <- ij <- aj*ai*p.j*p.k*p.i.il*(p.i.il-2*p.i)
    jl <- lj <- -aj*al*p.j*p.k*p.i.il*(p.l.il+2*p.l)
    kk <- ak^2*p.j*p.k.kil*((2*p.k-2+p.k.kil)*p.k*p.i.il+(1-p.k)*pi.kil-2*p.i.kil*p.k.kil)
    ki <- ik <- ak*ai*p.j*p.k.kil*(p.i.kil*(1-p.i)+2*p.i*p.k*p.i.il-2*p.i.kil^2+(p.i.kil-p.l.il)*p.k*p.i.il)
    kl <- lk <- ak*al*p.j*p.k.kil*(2*p.l*p.k*p.i.il-p.l*p.i.kil-2*p.l.kil*p.i.kil+p.k*p.i.il*(p.l.kil+p.l.il))
    ii <- 2*ai^2*p.j*p.k.kil*p.i.il*((p.i-1)*(p.i+p.i.kil)+p.i.kil^2+p.l.il*p.i.il)
    il <- li <- ai*al*p.j*p.k.kil*p.i.il*((p.i+p.i.kil-p.i.il)*(p.l+p.l.kil+p.l.il)+p.i.kil*p.l.kil+p.l*p.i-p.i.il*p.l.il)
    ll <- al^2*p.j*p.k.kil*p.i.il*((p.l+p.l.kil+p.l.il)^2-p.l.kil*(1-p.l.kil)-p.l*(1-p.l)-p.l.il*(1-p.l.il))
  }else if(identical(x,c(2,3,4,1))){
    jj <- aj^2*p.j*p.k*(1-2*p.j)*p.l.li
    jk <- kj <- aj*ak*p.j*p.k*(1-2*p.k)*p.l.li
    jl <- lj <- aj*al*p.j*p.k*p.l.li*(p.l.li-2*p.l)
    ji <- ij <- -aj*ai*p.j*p.k*p.l.li*(p.i.li+2*p.i)
    kk <- ak^2*p.j*p.k.kli*((2*p.k-2+p.k.kli)*p.k*p.l.li+(1-p.k)*pl.kli-2*p.l.kli*p.k.kli)
    kl <- lk <- ak*al*p.j*p.k.kli*(p.l.kli*(1-p.l)+2*p.l*p.k*p.l.li-2*p.l.kli^2+(p.l.kli-p.i.li)*p.k*p.l.li)
    ki <- ik <- ak*ai*p.j*p.k.kli*(2*p.i*p.k*p.l.li-p.i*p.l.kli-2*p.i.kli*p.l.kli+p.k*p.l.li*(p.i.kli+p.i.li))
    ll <- 2*al^2*p.j*p.k.kli*p.l.li*((p.l-1)*(p.l+p.l.kli)+p.l.kli^2+p.i.li*p.l.li)
    li <- il <- al*ai*p.j*p.k.kli*p.l.li*((p.l+p.l.kli-p.l.li)*(p.i+p.i.kli+p.i.li)+p.l.kli*p.i.kli+p.i*p.l-p.l.li*p.i.li)
    ii <- ai^2*p.j*p.k.kli*p.l.li*((p.i+p.i.kli+p.i.li)^2-p.i.kli*(1-p.i.kli)-p.i*(1-p.i)-p.i.li*(1-p.i.li))
  }else if(identical(x,c(2,4,1,3))){
    jj <- aj^2*p.j*p.l*(1-2*p.j)*p.i.ik
    jl <- lj <- aj*al*p.j*p.l*(1-2*p.l)*p.i.ik
    ji <- ij <- aj*ai*p.j*p.l*p.i.ik*(p.i.ik-2*p.i)
    jk <- kj <- -aj*ak*p.j*p.l*p.i.ik*(p.k.ik+2*p.k)
    ll <- al^2*p.j*p.l.lik*((2*p.l-2+p.l.lik)*p.l*p.i.ik+(1-p.l)*pi.lik-2*p.i.lik*p.l.lik)
    li <- il <- al*ai*p.j*p.l.lik*(p.i.lik*(1-p.i)+2*p.i*p.l*p.i.ik-2*p.i.lik^2+(p.i.lik-p.k.ik)*p.l*p.i.ik)
    lk <- kl <- al*ak*p.j*p.l.lik*(2*p.k*p.l*p.i.ik-p.k*p.i.lik-2*p.k.lik*p.i.lik+p.l*p.i.ik*(p.k.lik+p.k.ik))
    ii <- 2*ai^2*p.j*p.l.lik*p.i.ik*((p.i-1)*(p.i+p.i.lik)+p.i.lik^2+p.k.ik*p.i.ik)
    ik <- ki <- ai*ak*p.j*p.l.lik*p.i.ik*((p.i+p.i.lik-p.i.ik)*(p.k+p.k.lik+p.k.ik)+p.i.lik*p.k.lik+p.k*p.i-p.i.ik*p.k.ik)
    kk <- ak^2*p.j*p.l.lik*p.i.ik*((p.k+p.k.lik+p.k.ik)^2-p.k.lik*(1-p.k.lik)-p.k*(1-p.k)-p.k.ik*(1-p.k.ik))
  }else if(identical(x,c(2,4,3,1))){
    jj <- aj^2*p.j*p.l*(1-2*p.j)*p.k.ki
    jl <- lj <- aj*al*p.j*p.l*(1-2*p.l)*p.k.ki
    jk <- kj <- aj*ak*p.j*p.l*p.k.ki*(p.k.ki-2*p.k)
    ji <- ij <- -aj*ai*p.j*p.l*p.k.ki*(p.i.ki+2*p.i)
    ll <- al^2*p.j*p.l.lki*((2*p.l-2+p.l.lki)*p.l*p.k.ki+(1-p.l)*pk.lki-2*p.k.lki*p.l.lki)
    lk <- kl <- al*ak*p.j*p.l.lki*(p.k.lki*(1-p.k)+2*p.k*p.l*p.k.ki-2*p.k.lki^2+(p.k.lki-p.i.ki)*p.l*p.k.ki)
    li <- il <- al*ai*p.j*p.l.lki*(2*p.i*p.l*p.k.ki-p.i*p.k.lki-2*p.i.lki*p.k.lki+p.l*p.k.ki*(p.i.lki+p.i.ki))
    kk <- 2*ak^2*p.j*p.l.lki*p.k.ki*((p.k-1)*(p.k+p.k.lki)+p.k.lki^2+p.i.ki*p.k.ki)
    ki <- ik <- ak*ai*p.j*p.l.lki*p.k.ki*((p.k+p.k.lki-p.k.ki)*(p.i+p.i.lki+p.i.ki)+p.k.lki*p.i.lki+p.i*p.k-p.k.ki*p.i.ki)
    ii <- ai^2*p.j*p.l.lki*p.k.ki*((p.i+p.i.lki+p.i.ki)^2-p.i.lki*(1-p.i.lki)-p.i*(1-p.i)-p.i.ki*(1-p.i.ki))
  }else if(identical(x,c(3,1,2,4))){
    kk <- ak^2*p.k*p.i*(1-2*p.k)*p.j.jl
    ki <- ik <- ak*ai*p.k*p.i*(1-2*p.i)*p.j.jl
    kj <- jk <- ak*aj*p.k*p.i*p.j.jl*(p.j.jl-2*p.j)
    kl <- lk <- -ak*al*p.k*p.i*p.j.jl*(p.l.jl+2*p.l)
    ii <- ai^2*p.k*p.i.ijl*((2*p.i-2+p.i.ijl)*p.i*p.j.jl+(1-p.i)*pj.ijl-2*p.j.ijl*p.i.ijl)
    ij <- ji <- ai*aj*p.k*p.i.ijl*(p.j.ijl*(1-p.j)+2*p.j*p.i*p.j.jl-2*p.j.ijl^2+(p.j.ijl-p.l.jl)*p.i*p.j.jl)
    il <- li <- ai*al*p.k*p.i.ijl*(2*p.l*p.i*p.j.jl-p.l*p.j.ijl-2*p.l.ijl*p.j.ijl+p.i*p.j.jl*(p.l.ijl+p.l.jl))
    jj <- 2*aj^2*p.k*p.i.ijl*p.j.jl*((p.j-1)*(p.j+p.j.ijl)+p.j.ijl^2+p.l.jl*p.j.jl)
    jl <- lj <- aj*al*p.k*p.i.ijl*p.j.jl*((p.j+p.j.ijl-p.j.jl)*(p.l+p.l.ijl+p.l.jl)+p.j.ijl*p.l.ijl+p.l*p.j-p.j.jl*p.l.jl)
    ll <- al^2*p.k*p.i.ijl*p.j.jl*((p.l+p.l.ijl+p.l.jl)^2-p.l.ijl*(1-p.l.ijl)-p.l*(1-p.l)-p.l.jl*(1-p.l.jl))
  }else if(identical(x,c(3,1,4,2))){
    kk <- ak^2*p.k*p.i*(1-2*p.k)*p.l.lj
    ki <- ik <- ak*ai*p.k*p.i*(1-2*p.i)*p.l.lj
    kl <- lk <- ak*al*p.k*p.i*p.l.lj*(p.l.lj-2*p.l)
    kj <- jk <- -ak*aj*p.k*p.i*p.l.lj*(p.j.lj+2*p.j)
    ii <- ai^2*p.k*p.i.ilj*((2*p.i-2+p.i.ilj)*p.i*p.l.lj+(1-p.i)*pl.ilj-2*p.l.ilj*p.i.ilj)
    il <- li <- ai*al*p.k*p.i.ilj*(p.l.ilj*(1-p.l)+2*p.l*p.i*p.l.lj-2*p.l.ilj^2+(p.l.ilj-p.j.lj)*p.i*p.l.lj)
    ij <- ji <- ai*aj*p.k*p.i.ilj*(2*p.j*p.i*p.l.lj-p.j*p.l.ilj-2*p.j.ilj*p.l.ilj+p.i*p.l.lj*(p.j.ilj+p.j.lj))
    ll <- 2*al^2*p.k*p.i.ilj*p.l.lj*((p.l-1)*(p.l+p.l.ilj)+p.l.ilj^2+p.j.lj*p.l.lj)
    lj <- jl <- al*aj*p.k*p.i.ilj*p.l.lj*((p.l+p.l.ilj-p.l.lj)*(p.j+p.j.ilj+p.j.lj)+p.l.ilj*p.j.ilj+p.j*p.l-p.l.lj*p.j.lj)
    jj <- aj^2*p.k*p.i.ilj*p.l.lj*((p.j+p.j.ilj+p.j.lj)^2-p.j.ilj*(1-p.j.ilj)-p.j*(1-p.j)-p.j.lj*(1-p.j.lj))
  }else if(identical(x,c(3,2,1,4))){
    kk <- ak^2*p.k*p.j*(1-2*p.k)*p.i.il
    kj <- jk <- ak*aj*p.k*p.j*(1-2*p.j)*p.i.il
    ki <- ik <- ak*ai*p.k*p.j*p.i.il*(p.i.il-2*p.i)
    kl <- lk <- -ak*al*p.k*p.j*p.i.il*(p.l.il+2*p.l)
    jj <- aj^2*p.k*p.j.jil*((2*p.j-2+p.j.jil)*p.j*p.i.il+(1-p.j)*pi.jil-2*p.i.jil*p.j.jil)
    ji <- ij <- aj*ai*p.k*p.j.jil*(p.i.jil*(1-p.i)+2*p.i*p.j*p.i.il-2*p.i.jil^2+(p.i.jil-p.l.il)*p.j*p.i.il)
    jl <- lj <- aj*al*p.k*p.j.jil*(2*p.l*p.j*p.i.il-p.l*p.i.jil-2*p.l.jil*p.i.jil+p.j*p.i.il*(p.l.jil+p.l.il))
    ii <- 2*ai^2*p.k*p.j.jil*p.i.il*((p.i-1)*(p.i+p.i.jil)+p.i.jil^2+p.l.il*p.i.il)
    il <- li <- ai*al*p.k*p.j.jil*p.i.il*((p.i+p.i.jil-p.i.il)*(p.l+p.l.jil+p.l.il)+p.i.jil*p.l.jil+p.l*p.i-p.i.il*p.l.il)
    ll <- al^2*p.k*p.j.jil*p.i.il*((p.l+p.l.jil+p.l.il)^2-p.l.jil*(1-p.l.jil)-p.l*(1-p.l)-p.l.il*(1-p.l.il))
  }else if(identical(x,c(3,2,4,1))){
    kk <- ak^2*p.k*p.j*(1-2*p.k)*p.l.li
    kj <- jk <- ak*aj*p.k*p.j*(1-2*p.j)*p.l.li
    kl <- lk <- ak*al*p.k*p.j*p.l.li*(p.l.li-2*p.l)
    ki <- ik <- -ak*ai*p.k*p.j*p.l.li*(p.i.li+2*p.i)
    jj <- aj^2*p.k*p.j.jli*((2*p.j-2+p.j.jli)*p.j*p.l.li+(1-p.j)*pl.jli-2*p.l.jli*p.j.jli)
    jl <- lj <- aj*al*p.k*p.j.jli*(p.l.jli*(1-p.l)+2*p.l*p.j*p.l.li-2*p.l.jli^2+(p.l.jli-p.i.li)*p.j*p.l.li)
    ji <- ij <- aj*ai*p.k*p.j.jli*(2*p.i*p.j*p.l.li-p.i*p.l.jli-2*p.i.jli*p.l.jli+p.j*p.l.li*(p.i.jli+p.i.li))
    ll <- 2*al^2*p.k*p.j.jli*p.l.li*((p.l-1)*(p.l+p.l.jli)+p.l.jli^2+p.i.li*p.l.li)
    li <- il <- al*ai*p.k*p.j.jli*p.l.li*((p.l+p.l.jli-p.l.li)*(p.i+p.i.jli+p.i.li)+p.l.jli*p.i.jli+p.i*p.l-p.l.li*p.i.li)
    ii <- ai^2*p.k*p.j.jli*p.l.li*((p.i+p.i.jli+p.i.li)^2-p.i.jli*(1-p.i.jli)-p.i*(1-p.i)-p.i.li*(1-p.i.li))
  }else if(identical(x,c(3,4,1,2))){
    kk <- ak^2*p.k*p.l*(1-2*p.k)*p.i.ij
    kl <- lk <- ak*al*p.k*p.l*(1-2*p.l)*p.i.ij
    ki <- ik <- ak*ai*p.k*p.l*p.i.ij*(p.i.ij-2*p.i)
    kj <- jk <- -ak*aj*p.k*p.l*p.i.ij*(p.j.ij+2*p.j)
    ll <- al^2*p.k*p.l.lij*((2*p.l-2+p.l.lij)*p.l*p.i.ij+(1-p.l)*pi.lij-2*p.i.lij*p.l.lij)
    li <- il <- al*ai*p.k*p.l.lij*(p.i.lij*(1-p.i)+2*p.i*p.l*p.i.ij-2*p.i.lij^2+(p.i.lij-p.j.ij)*p.l*p.i.ij)
    lj <- jl <- al*aj*p.k*p.l.lij*(2*p.j*p.l*p.i.ij-p.j*p.i.lij-2*p.j.lij*p.i.lij+p.l*p.i.ij*(p.j.lij+p.j.ij))
    ii <- 2*ai^2*p.k*p.l.lij*p.i.ij*((p.i-1)*(p.i+p.i.lij)+p.i.lij^2+p.j.ij*p.i.ij)
    ij <- ji <- ai*aj*p.k*p.l.lij*p.i.ij*((p.i+p.i.lij-p.i.ij)*(p.j+p.j.lij+p.j.ij)+p.i.lij*p.j.lij+p.j*p.i-p.i.ij*p.j.ij)
    jj <- aj^2*p.k*p.l.lij*p.i.ij*((p.j+p.j.lij+p.j.ij)^2-p.j.lij*(1-p.j.lij)-p.j*(1-p.j)-p.j.ij*(1-p.j.ij))
  }else if(identical(x,c(3,4,2,1))){
    kk <- ak^2*p.k*p.l*(1-2*p.k)*p.j.ji
    kl <- lk <- ak*al*p.k*p.l*(1-2*p.l)*p.j.ji
    kj <- jk <- ak*aj*p.k*p.l*p.j.ji*(p.j.ji-2*p.j)
    ki <- ik <- -ak*ai*p.k*p.l*p.j.ji*(p.i.ji+2*p.i)
    ll <- al^2*p.k*p.l.lji*((2*p.l-2+p.l.lji)*p.l*p.j.ji+(1-p.l)*pj.lji-2*p.j.lji*p.l.lji)
    lj <- jl <- al*aj*p.k*p.l.lji*(p.j.lji*(1-p.j)+2*p.j*p.l*p.j.ji-2*p.j.lji^2+(p.j.lji-p.i.ji)*p.l*p.j.ji)
    li <- il <- al*ai*p.k*p.l.lji*(2*p.i*p.l*p.j.ji-p.i*p.j.lji-2*p.i.lji*p.j.lji+p.l*p.j.ji*(p.i.lji+p.i.ji))
    jj <- 2*aj^2*p.k*p.l.lji*p.j.ji*((p.j-1)*(p.j+p.j.lji)+p.j.lji^2+p.i.ji*p.j.ji)
    ji <- ij <- aj*ai*p.k*p.l.lji*p.j.ji*((p.j+p.j.lji-p.j.ji)*(p.i+p.i.lji+p.i.ji)+p.j.lji*p.i.lji+p.i*p.j-p.j.ji*p.i.ji)
    ii <- ai^2*p.k*p.l.lji*p.j.ji*((p.i+p.i.lji+p.i.ji)^2-p.i.lji*(1-p.i.lji)-p.i*(1-p.i)-p.i.ji*(1-p.i.ji))
  }else if(identical(x,c(4,1,2,3))){
    ll <- al^2*p.l*p.i*(1-2*p.l)*p.j.jk
    li <- il <- al*ai*p.l*p.i*(1-2*p.i)*p.j.jk
    lj <- jl <- al*aj*p.l*p.i*p.j.jk*(p.j.jk-2*p.j)
    lk <- kl <- -al*ak*p.l*p.i*p.j.jk*(p.k.jk+2*p.k)
    ii <- ai^2*p.l*p.i.ijk*((2*p.i-2+p.i.ijk)*p.i*p.j.jk+(1-p.i)*pj.ijk-2*p.j.ijk*p.i.ijk)
    ij <- ji <- ai*aj*p.l*p.i.ijk*(p.j.ijk*(1-p.j)+2*p.j*p.i*p.j.jk-2*p.j.ijk^2+(p.j.ijk-p.k.jk)*p.i*p.j.jk)
    ik <- ki <- ai*ak*p.l*p.i.ijk*(2*p.k*p.i*p.j.jk-p.k*p.j.ijk-2*p.k.ijk*p.j.ijk+p.i*p.j.jk*(p.k.ijk+p.k.jk))
    jj <- 2*aj^2*p.l*p.i.ijk*p.j.jk*((p.j-1)*(p.j+p.j.ijk)+p.j.ijk^2+p.k.jk*p.j.jk)
    jk <- kj <- aj*ak*p.l*p.i.ijk*p.j.jk*((p.j+p.j.ijk-p.j.jk)*(p.k+p.k.ijk+p.k.jk)+p.j.ijk*p.k.ijk+p.k*p.j-p.j.jk*p.k.jk)
    kk <- ak^2*p.l*p.i.ijk*p.j.jk*((p.k+p.k.ijk+p.k.jk)^2-p.k.ijk*(1-p.k.ijk)-p.k*(1-p.k)-p.k.jk*(1-p.k.jk))
  }else if(identical(x,c(4,1,3,2))){
    ll <- al^2*p.l*p.i*(1-2*p.l)*p.k.kj
    li <- il <- al*ai*p.l*p.i*(1-2*p.i)*p.k.kj
    lk <- kl <- al*ak*p.l*p.i*p.k.kj*(p.k.kj-2*p.k)
    lj <- jl <- -al*aj*p.l*p.i*p.k.kj*(p.j.kj+2*p.j)
    ii <- ai^2*p.l*p.i.ikj*((2*p.i-2+p.i.ikj)*p.i*p.k.kj+(1-p.i)*pk.ikj-2*p.k.ikj*p.i.ikj)
    ik <- ki <- ai*ak*p.l*p.i.ikj*(p.k.ikj*(1-p.k)+2*p.k*p.i*p.k.kj-2*p.k.ikj^2+(p.k.ikj-p.j.kj)*p.i*p.k.kj)
    ij <- ji <- ai*aj*p.l*p.i.ikj*(2*p.j*p.i*p.k.kj-p.j*p.k.ikj-2*p.j.ikj*p.k.ikj+p.i*p.k.kj*(p.j.ikj+p.j.kj))
    kk <- 2*ak^2*p.l*p.i.ikj*p.k.kj*((p.k-1)*(p.k+p.k.ikj)+p.k.ikj^2+p.j.kj*p.k.kj)
    kj <- jk <- ak*aj*p.l*p.i.ikj*p.k.kj*((p.k+p.k.ikj-p.k.kj)*(p.j+p.j.ikj+p.j.kj)+p.k.ikj*p.j.ikj+p.j*p.k-p.k.kj*p.j.kj)
    jj <- aj^2*p.l*p.i.ikj*p.k.kj*((p.j+p.j.ikj+p.j.kj)^2-p.j.ikj*(1-p.j.ikj)-p.j*(1-p.j)-p.j.kj*(1-p.j.kj))
  }else if(identical(x,c(4,2,1,3))){
    ll <- al^2*p.l*p.j*(1-2*p.l)*p.i.ik
    lj <- jl <- al*aj*p.l*p.j*(1-2*p.j)*p.i.ik
    li <- il <- al*ai*p.l*p.j*p.i.ik*(p.i.ik-2*p.i)
    lk <- kl <- -al*ak*p.l*p.j*p.i.ik*(p.k.ik+2*p.k)
    jj <- aj^2*p.l*p.j.jik*((2*p.j-2+p.j.jik)*p.j*p.i.ik+(1-p.j)*pi.jik-2*p.i.jik*p.j.jik)
    ji <- ij <- aj*ai*p.l*p.j.jik*(p.i.jik*(1-p.i)+2*p.i*p.j*p.i.ik-2*p.i.jik^2+(p.i.jik-p.k.ik)*p.j*p.i.ik)
    jk <- kj <- aj*ak*p.l*p.j.jik*(2*p.k*p.j*p.i.ik-p.k*p.i.jik-2*p.k.jik*p.i.jik+p.j*p.i.ik*(p.k.jik+p.k.ik))
    ii <- 2*ai^2*p.l*p.j.jik*p.i.ik*((p.i-1)*(p.i+p.i.jik)+p.i.jik^2+p.k.ik*p.i.ik)
    ik <- ki <- ai*ak*p.l*p.j.jik*p.i.ik*((p.i+p.i.jik-p.i.ik)*(p.k+p.k.jik+p.k.ik)+p.i.jik*p.k.jik+p.k*p.i-p.i.ik*p.k.ik)
    kk <- ak^2*p.l*p.j.jik*p.i.ik*((p.k+p.k.jik+p.k.ik)^2-p.k.jik*(1-p.k.jik)-p.k*(1-p.k)-p.k.ik*(1-p.k.ik))
  }else if(identical(x,c(4,2,3,1))){
    ll <- al^2*p.l*p.j*(1-2*p.l)*p.k.ki
    lj <- jl <- al*aj*p.l*p.j*(1-2*p.j)*p.k.ki
    lk <- kl <- al*ak*p.l*p.j*p.k.ki*(p.k.ki-2*p.k)
    li <- il <- -al*ai*p.l*p.j*p.k.ki*(p.i.ki+2*p.i)
    jj <- aj^2*p.l*p.j.jki*((2*p.j-2+p.j.jki)*p.j*p.k.ki+(1-p.j)*pk.jki-2*p.k.jki*p.j.jki)
    jk <- kj <- aj*ak*p.l*p.j.jki*(p.k.jki*(1-p.k)+2*p.k*p.j*p.k.ki-2*p.k.jki^2+(p.k.jki-p.i.ki)*p.j*p.k.ki)
    ji <- ij <- aj*ai*p.l*p.j.jki*(2*p.i*p.j*p.k.ki-p.i*p.k.jki-2*p.i.jki*p.k.jki+p.j*p.k.ki*(p.i.jki+p.i.ki))
    kk <- 2*ak^2*p.l*p.j.jki*p.k.ki*((p.k-1)*(p.k+p.k.jki)+p.k.jki^2+p.i.ki*p.k.ki)
    ki <- ik <- ak*ai*p.l*p.j.jki*p.k.ki*((p.k+p.k.jki-p.k.ki)*(p.i+p.i.jki+p.i.ki)+p.k.jki*p.i.jki+p.i*p.k-p.k.ki*p.i.ki)
    ii <- ai^2*p.l*p.j.jki*p.k.ki*((p.i+p.i.jki+p.i.ki)^2-p.i.jki*(1-p.i.jki)-p.i*(1-p.i)-p.i.ki*(1-p.i.ki))
  }else if(identical(x,c(4,3,1,2))){
    ll <- al^2*p.l*p.k*(1-2*p.l)*p.i.ij
    lk <- kl <- al*ak*p.l*p.k*(1-2*p.k)*p.i.ij
    li <- il <- al*ai*p.l*p.k*p.i.ij*(p.i.ij-2*p.i)
    lj <- jl <- -al*aj*p.l*p.k*p.i.ij*(p.j.ij+2*p.j)
    kk <- ak^2*p.l*p.k.kij*((2*p.k-2+p.k.kij)*p.k*p.i.ij+(1-p.k)*pi.kij-2*p.i.kij*p.k.kij)
    ki <- ik <- ak*ai*p.l*p.k.kij*(p.i.kij*(1-p.i)+2*p.i*p.k*p.i.ij-2*p.i.kij^2+(p.i.kij-p.j.ij)*p.k*p.i.ij)
    kj <- jk <- ak*aj*p.l*p.k.kij*(2*p.j*p.k*p.i.ij-p.j*p.i.kij-2*p.j.kij*p.i.kij+p.k*p.i.ij*(p.j.kij+p.j.ij))
    ii <- 2*ai^2*p.l*p.k.kij*p.i.ij*((p.i-1)*(p.i+p.i.kij)+p.i.kij^2+p.j.ij*p.i.ij)
    ij <- ji <- ai*aj*p.l*p.k.kij*p.i.ij*((p.i+p.i.kij-p.i.ij)*(p.j+p.j.kij+p.j.ij)+p.i.kij*p.j.kij+p.j*p.i-p.i.ij*p.j.ij)
    jj <- aj^2*p.l*p.k.kij*p.i.ij*((p.j+p.j.kij+p.j.ij)^2-p.j.kij*(1-p.j.kij)-p.j*(1-p.j)-p.j.ij*(1-p.j.ij))
  }else if(identical(x,c(4,3,2,1))){
    ll <- al^2*p.l*p.k*(1-2*p.l)*p.j.ji
    lk <- kl <- al*ak*p.l*p.k*(1-2*p.k)*p.j.ji
    lj <- jl <- al*aj*p.l*p.k*p.j.ji*(p.j.ji-2*p.j)
    li <- il <- -al*ai*p.l*p.k*p.j.ji*(p.i.ji+2*p.i)
    kk <- ak^2*p.l*p.k.kji*((2*p.k-2+p.k.kji)*p.k*p.j.ji+(1-p.k)*pj.kji-2*p.j.kji*p.k.kji)
    kj <- jk <- ak*aj*p.l*p.k.kji*(p.j.kji*(1-p.j)+2*p.j*p.k*p.j.ji-2*p.j.kji^2+(p.j.kji-p.i.ji)*p.k*p.j.ji)
    ki <- ik <- ak*ai*p.l*p.k.kji*(2*p.i*p.k*p.j.ji-p.i*p.j.kji-2*p.i.kji*p.j.kji+p.k*p.j.ji*(p.i.kji+p.i.ji))
    jj <- 2*aj^2*p.l*p.k.kji*p.j.ji*((p.j-1)*(p.j+p.j.kji)+p.j.kji^2+p.i.ji*p.j.ji)
    ji <- ij <- aj*ai*p.l*p.k.kji*p.j.ji*((p.j+p.j.kji-p.j.ji)*(p.i+p.i.kji+p.i.ji)+p.j.kji*p.i.kji+p.i*p.j-p.j.ji*p.i.ji)
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


Ib_theta_2PL_4_mole <- function(theta,a,d){
  all.x <- list(c(1,2,3,4),c(1,2,4,3),c(1,3,2,4),c(1,3,4,2),c(1,4,2,3),c(1,4,3,2),
                c(2,1,3,4),c(2,1,4,3),c(2,3,1,4),c(2,3,4,1),c(2,4,1,3),c(2,4,3,1),
                c(3,1,2,4),c(3,1,4,2),c(3,2,1,4),c(3,2,4,1),c(3,4,1,2),c(3,4,2,1),
                c(4,1,2,3),c(4,1,3,2),c(4,2,1,3),c(4,2,3,1),c(4,3,1,2),c(4,3,2,1))
  duiying <- matrix(c(1,3,
                      2,5,
                      4,6,
                      7,9,
                      8,11,
                      10,12,
                      13,15,
                      14,17,
                      16,18,
                      19,21,
                      20,23,
                      22,24),ncol = 2,byrow = TRUE)
  ret <- 0
  for(x in 1:length(duiying)){
    g1 <- grad_theta_2PL_4_mole(theta,a,d,x=all.x[[duiying[x,1]]])
    mid1 <- tcrossprod(g1,g1)/px_2PL_4_mole(theta,a,d,x=all.x[[duiying[x,1]]])-hes_theta_2PL_4_mole(theta,a,d,x=all.x[[duiying[x,1]]])
    g2 <- grad_theta_2PL_4_mole(theta,a,d,x=all.x[[duiying[x,2]]])
    mid2 <- tcrossprod(g2,g2)/px_2PL_4_mole(theta,a,d,x=all.x[[duiying[x,2]]])-hes_theta_2PL_4_mole(theta,a,d,x=all.x[[duiying[x,2]]])
    ret <- ret + mid1+mid2
  }
  ret
}


test.info_theta_2PL_4_mole <- function(theta,item.par,BID){
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
    complete.info.b[dim.b,dim.b] <- info[[b]] <- Ib_theta_2PL_4_mole(theta.b,a.b,d.b)
    I <- I+complete.info.b
    comp.info[[b]] <- complete.info.b
  }
  se <- sqrt(diag(solve(Reduce("+",comp.info))))
  list(comp.info=comp.info,se =se, info.b=info,I=I)
}
