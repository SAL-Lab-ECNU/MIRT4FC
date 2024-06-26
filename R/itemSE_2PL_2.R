###################### SE计算 ################################

# P_b(x|theta)
## 入参：px以block为单位计算
## theta=作答者能力值，a=题目的a参数，d=（ a * b），xx=作答模式，会在Ib中遍历6种作答模式
px <- function(theta,a,d,xx=c(1,2,3)){ #a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)
  z <- exp(a*theta-d)
  prod(z[xx[-3]])/(sum(z)*sum(z[xx[-1]]))
}

grad_item <- function(theta,a,d,x=c(1,2,3)){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)
  z <- exp(a*theta-d)          #1*3
  thetai <- theta[1]
  thetaj <- theta[2]
  thetak <- theta[3]
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
  if(identical(x,c(1,2,3))){       #顺序为1,2,3时的一阶导，下方同理
    g.ia <- thetai*p.i*p.j
    g.ja <- thetaj*p.i*p.j.jk*(p.k.jk-p.j)
    g.ka <- -thetak*p.i*p.j.jk*(p.k+p.k.jk)
    g.id <- -p.i*p.j
    g.jd <- -p.i*p.j.jk*(p.k.jk-p.j)
    #g.kd <- p.i*p.j.jk*(p.k+p.k.jk)
  }else if(identical(x,c(2,3,1))){
    g.ja <- thetaj*p.j*p.k
    g.ka <- thetak*p.j*p.k.ki*(p.i.ki-p.k)
    g.ia <- -thetai*p.j*p.k.ki*(p.i+p.i.ki)
    g.jd <- -p.j*p.k
    #g.kd <- -p.j*p.k.ki*(p.i.ki-p.k)
    g.id <- p.j*p.k.ki*(p.i+p.i.ki)
  }else if(identical(x,c(1,3,2))){
    g.ia <- thetai*p.i*p.k
    g.ka <- thetak*p.i*p.k.kj*(p.j.kj-p.k)
    g.ja <- -thetaj*p.i*p.k.kj*(p.j+p.j.kj)
    g.id <- -p.i*p.k
    #g.kd <- -p.i*p.k.kj*(p.j.kj-p.k)
    g.jd <- p.i*p.k.kj*(p.j+p.j.kj)
  }else if(identical(x,c(2,1,3))){
    g.ja <- thetaj*p.j*p.i
    g.ia <- thetai*p.j*p.i.ik*(p.k.ik-p.i)
    g.ka <- -thetak*p.j*p.i.ik*(p.k+p.k.ik)
    g.jd <- -p.j*p.i
    g.id <- -p.j*p.i.ik*(p.k.ik-p.i)
    #g.kd <- p.j*p.i.ik*(p.k+p.k.ik)
  }else if(identical(x,c(3,1,2))){
    g.ka <- thetak*p.k*p.i
    g.ia <- thetai*p.k*p.i.ij*(p.j.ij-p.i)
    g.ja <- -thetaj*p.k*p.i.ij*(p.j+p.j.ij)
    #g.kd <- -p.k*p.i
    g.id <- -p.k*p.i.ij*(p.j.ij-p.i)
    g.jd <- p.k*p.i.ij*(p.j+p.j.ij)
  }else if(identical(x,c(3,2,1))){
    g.ka <- thetak*p.k*p.j
    g.ja <- thetaj*p.k*p.j.ji*(p.i.ji-p.j)
    g.ia <- -thetai*p.k*p.j.ji*(p.i+p.i.ji)
    #g.kd <- -p.k*p.j
    g.jd <- -p.k*p.j.ji*(p.i.ji-p.j)
    g.id <- p.k*p.j.ji*(p.i+p.i.ji)
  }
  
  c(g.ia,g.ja,g.ka,g.id,g.jd)
}

hes_item <- function(theta,a,d,x=c(1,2,3)){#一个block内的二阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)
  h <- matrix(NA,5,5)                 #初始化海森矩阵，ad两参数各3个，因此为6*6
  z <- exp(a*theta-d)
  thetai <- theta[1]
  thetaj <- theta[2]
  thetak <- theta[3]
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
  
  if(identical(x,c(1,2,3))){                #1,2,3顺序下的海森矩阵详细计算公式，下同
    iaia <- thetai^2*p.i*p.j*(1-2*p.i)
    jaja <- -thetaj^2*p.i*p.j*p.j.jk*(1-2*p.j+p.k.jk)-thetaj^2*p.i*p.j.jk*p.k.jk*(p.j-p.k.jk+p.j.jk)
    kaka <- thetak^2*p.i*p.j.jk*((p.k+p.k.jk)^2+p.k^2-p.k-p.j.jk*p.k.jk)
    idid <- p.i*p.j*(1-2*p.i)
    jdjd <- -p.i*p.j*p.j.jk*(1-2*p.j+p.k.jk)-p.i*p.j.jk*p.k.jk*(p.j-p.k.jk+p.j.jk)
    #kdkd <- p.i*p.j.jk*((p.k+p.k.jk)^2+p.k^2-p.k-p.j.jk*p.k.jk)
    iaja <- jaia <- thetai*thetaj*p.i*p.j*(1-2*p.j)
    iaka <- kaia <- -2*thetai*thetak*p.i*p.j*p.k
    iaid <- idia <- -thetai*p.i*p.j*(1-2*p.i)
    iajd <- jdia <- -thetai*p.i*p.j*(1-2*p.j)
    #iakd <- kdia <- 2*thetai*p.i*p.j*p.k
    jaka <- kaja <- thetaj*thetak*p.i*p.j.jk*(-1*p.k*p.k.jk+2*p.k*p.j-p.k.jk*p.k.jk+p.k.jk*p.j+p.j.jk*p.k.jk)
    jaid <- idja <- -thetaj*p.i*p.j*(1-2*p.j)
    jajd <- jdja <- thetaj*p.i*p.j*p.j.jk*(1-2*p.j+p.k.jk)+thetaj*p.i*p.j.jk*p.k.jk*(p.j-p.k.jk+p.j.jk)
    #jakd <- kdja <- -thetaj*p.i*p.j.jk*(-1*p.k*p.k.jk+2*p.k*p.j-p.k.jk*p.k.jk+p.k.jk*p.j+p.j.jk*p.k.jk)
    kaid <- idka <- 2*thetak*p.i*p.j*p.k
    kajd <- jdka <- -thetak*p.i*p.j.jk*(-1*p.k*p.k.jk+2*p.k*p.j-p.k.jk*p.k.jk+p.k.jk*p.j+p.j.jk*p.k.jk)
    #kakd <- kdka <- -thetak*p.i*p.j.jk*((p.k+p.k.jk)^2+p.k^2-p.k-p.j.jk*p.k.jk)
    idjd <- jdid <- p.i*p.j*(1-2*p.j)
    #idkd <- kdid <- -2*p.i*p.j*p.k
    #jdkd <- kdjd <- p.i*p.j.jk*(-1*p.k*p.k.jk+2*p.k*p.j-p.k.jk*p.k.jk+p.k.jk*p.j+p.j.jk*p.k.jk)
  }else if(identical(x,c(2,3,1))){
    jaja <- thetaj^2*p.j*p.k*(1-2*p.j)
    kaka <- -thetak^2*p.j*p.k*p.k.ki*(1-2*p.k+p.i.ki)-thetak^2*p.j*p.k.ki*p.i.ki*(p.k-p.i.ki+p.k.ki)
    iaia <- thetai^2*p.j*p.k.ki*((p.i+p.i.ki)^2+p.i^2-p.i-p.k.ki*p.i.ki)
    jdjd <- p.j*p.k*(1-2*p.j)
    #kdkd <- -p.j*p.k*p.k.ki*(1-2*p.k+p.i.ki)-p.j*p.k.ki*p.i.ki*(p.k-p.i.ki+p.k.ki)
    idid <- p.j*p.k.ki*((p.i+p.i.ki)^2+p.i^2-p.i-p.k.ki*p.i.ki)
    jaka <- kaja <- thetaj*thetak*p.j*p.k*(1-2*p.k)
    jaia <- iaja <- -2*thetaj*thetai*p.j*p.k*p.i
    jajd <- jdja <- -thetaj*p.j*p.k*(1-2*p.j)
    #jakd <- kdja <- -thetaj*p.j*p.k*(1-2*p.k)
    jaid <- idja <- 2*thetaj*p.j*p.k*p.i
    kaia <- iaka <- thetak*thetai*p.j*p.k.ki*(-1*p.i*p.i.ki+2*p.i*p.k-p.i.ki*p.i.ki+p.i.ki*p.k+p.k.ki*p.i.ki)
    kajd <- jdka <- -thetak*p.j*p.k*(1-2*p.k)
    #kakd <- kdka <- thetak*p.j*p.k*p.k.ki*(1-2*p.k+p.i.ki)+thetak*p.j*p.k.ki*p.i.ki*(p.k-p.i.ki+p.k.ki)
    kaid <- idka <- -thetak*p.j*p.k.ki*(-1*p.i*p.i.ki+2*p.i*p.k-p.i.ki*p.i.ki+p.i.ki*p.k+p.k.ki*p.i.ki)
    iajd <- jdia <- 2*thetai*p.j*p.k*p.i
    #iakd <- kdia <- -thetai*p.j*p.k.ki*(-1*p.i*p.i.ki+2*p.i*p.k-p.i.ki*p.i.ki+p.i.ki*p.k+p.k.ki*p.i.ki)
    iaid <- idia <- -thetai*p.j*p.k.ki*((p.i+p.i.ki)^2+p.i^2-p.i-p.k.ki*p.i.ki)
    #jdkd <- kdjd <- p.j*p.k*(1-2*p.k)
    jdid <- idjd <- -2*p.j*p.k*p.i
    #kdid <- idkd <- p.j*p.k.ki*(-1*p.i*p.i.ki+2*p.i*p.k-p.i.ki*p.i.ki+p.i.ki*p.k+p.k.ki*p.i.ki)
  }else if(identical(x,c(1,3,2))){
    iaia <- thetai^2*p.i*p.k*(1-2*p.i)
    kaka <- -thetak^2*p.i*p.k*p.k.kj*(1-2*p.k+p.j.kj)-thetak^2*p.i*p.k.kj*p.j.kj*(p.k-p.j.kj+p.k.kj)
    jaja <- thetaj^2*p.i*p.k.kj*((p.j+p.j.kj)^2+p.j^2-p.j-p.k.kj*p.j.kj)
    idid <- p.i*p.k*(1-2*p.i)
    #kdkd <- -p.i*p.k*p.k.kj*(1-2*p.k+p.j.kj)-p.i*p.k.kj*p.j.kj*(p.k-p.j.kj+p.k.kj)
    jdjd <- p.i*p.k.kj*((p.j+p.j.kj)^2+p.j^2-p.j-p.k.kj*p.j.kj)
    iaka <- kaia <- thetai*thetak*p.i*p.k*(1-2*p.k)
    iaja <- jaia <- -2*thetai*thetaj*p.i*p.k*p.j
    iaid <- idia <- -thetai*p.i*p.k*(1-2*p.i)
    #iakd <- kdia <- -thetai*p.i*p.k*(1-2*p.k)
    iajd <- jdia <- 2*thetai*p.i*p.k*p.j
    kaja <- jaka <- thetak*thetaj*p.i*p.k.kj*(-1*p.j*p.j.kj+2*p.j*p.k-p.j.kj*p.j.kj+p.j.kj*p.k+p.k.kj*p.j.kj)
    kaid <- idka <- -thetak*p.i*p.k*(1-2*p.k)
    #kakd <- kdka <- thetak*p.i*p.k*p.k.kj*(1-2*p.k+p.j.kj)+thetak*p.i*p.k.kj*p.j.kj*(p.k-p.j.kj+p.k.kj)
    kajd <- jdka <- -thetak*p.i*p.k.kj*(-1*p.j*p.j.kj+2*p.j*p.k-p.j.kj*p.j.kj+p.j.kj*p.k+p.k.kj*p.j.kj)
    jaid <- idja <- 2*thetaj*p.i*p.k*p.j
    #jakd <- kdja <- -thetaj*p.i*p.k.kj*(-1*p.j*p.j.kj+2*p.j*p.k-p.j.kj*p.j.kj+p.j.kj*p.k+p.k.kj*p.j.kj)
    jajd <- jdja <- -thetaj*p.i*p.k.kj*((p.j+p.j.kj)^2+p.j^2-p.j-p.k.kj*p.j.kj)
    #idkd <- kdid <- p.i*p.k*(1-2*p.k)
    idjd <- jdid <- -2*p.i*p.k*p.j
    #kdjd <- jdkd <- p.i*p.k.kj*(-1*p.j*p.j.kj+2*p.j*p.k-p.j.kj*p.j.kj+p.j.kj*p.k+p.k.kj*p.j.kj)
  }else if(identical(x,c(2,1,3))){
    jaja <- thetaj^2*p.j*p.i*(1-2*p.j)
    iaia <- -thetai^2*p.j*p.i*p.i.ik*(1-2*p.i+p.k.ik)-thetai^2*p.j*p.i.ik*p.k.ik*(p.i-p.k.ik+p.i.ik)
    kaka <- thetak^2*p.j*p.i.ik*((p.k+p.k.ik)^2+p.k^2-p.k-p.i.ik*p.k.ik)
    jdjd <- p.j*p.i*(1-2*p.j)
    idid <- -p.j*p.i*p.i.ik*(1-2*p.i+p.k.ik)-p.j*p.i.ik*p.k.ik*(p.i-p.k.ik+p.i.ik)
    #kdkd <- p.j*p.i.ik*((p.k+p.k.ik)^2+p.k^2-p.k-p.i.ik*p.k.ik)
    jaia <- iaja <- thetaj*thetai*p.j*p.i*(1-2*p.i)
    jaka <- kaja <- -2*thetaj*thetak*p.j*p.i*p.k
    jajd <- jdja <- -thetaj*p.j*p.i*(1-2*p.j)
    jaid <- idja <- -thetaj*p.j*p.i*(1-2*p.i)
    #jakd <- kdja <- 2*thetaj*p.j*p.i*p.k
    iaka <- kaia <- thetai*thetak*p.j*p.i.ik*(-1*p.k*p.k.ik+2*p.k*p.i-p.k.ik*p.k.ik+p.k.ik*p.i+p.i.ik*p.k.ik)
    iajd <- jdia <- -thetai*p.j*p.i*(1-2*p.i)
    iaid <- idia <- thetai*p.j*p.i*p.i.ik*(1-2*p.i+p.k.ik)+thetai*p.j*p.i.ik*p.k.ik*(p.i-p.k.ik+p.i.ik)
    #iakd <- kdia <- -thetai*p.j*p.i.ik*(-1*p.k*p.k.ik+2*p.k*p.i-p.k.ik*p.k.ik+p.k.ik*p.i+p.i.ik*p.k.ik)
    kajd <- jdka <- 2*thetak*p.j*p.i*p.k
    kaid <- idka <- -thetak*p.j*p.i.ik*(-1*p.k*p.k.ik+2*p.k*p.i-p.k.ik*p.k.ik+p.k.ik*p.i+p.i.ik*p.k.ik)
    #kakd <- kdka <- -thetak*p.j*p.i.ik*((p.k+p.k.ik)^2+p.k^2-p.k-p.i.ik*p.k.ik)
    jdid <- idjd <- p.j*p.i*(1-2*p.i)
    #jdkd <- kdjd <- -2*p.j*p.i*p.k
    #idkd <- kdid <- p.j*p.i.ik*(-1*p.k*p.k.ik+2*p.k*p.i-p.k.ik*p.k.ik+p.k.ik*p.i+p.i.ik*p.k.ik)
  }else if(identical(x,c(3,1,2))){
    kaka <- thetak^2*p.k*p.i*(1-2*p.k)
    iaia <- -thetai^2*p.k*p.i*p.i.ij*(1-2*p.i+p.j.ij)-thetai^2*p.k*p.i.ij*p.j.ij*(p.i-p.j.ij+p.i.ij)
    jaja <- thetaj^2*p.k*p.i.ij*((p.j+p.j.ij)^2+p.j^2-p.j-p.i.ij*p.j.ij)
    #kdkd <- p.k*p.i*(1-2*p.k)
    idid <- -p.k*p.i*p.i.ij*(1-2*p.i+p.j.ij)-p.k*p.i.ij*p.j.ij*(p.i-p.j.ij+p.i.ij)
    jdjd <- p.k*p.i.ij*((p.j+p.j.ij)^2+p.j^2-p.j-p.i.ij*p.j.ij)
    kaia <- iaka <- thetak*thetai*p.k*p.i*(1-2*p.i)
    kaja <- jaka <- -2*thetak*thetaj*p.k*p.i*p.j
    #kakd <- kdka <- -thetak*p.k*p.i*(1-2*p.k)
    kaid <- idka <- -thetak*p.k*p.i*(1-2*p.i)
    kajd <- jdka <- 2*thetak*p.k*p.i*p.j
    iaja <- jaia <- thetai*thetaj*p.k*p.i.ij*(-1*p.j*p.j.ij+2*p.j*p.i-p.j.ij*p.j.ij+p.j.ij*p.i+p.i.ij*p.j.ij)
    #iakd <- kdia <- -thetai*p.k*p.i*(1-2*p.i)
    iaid <- idia <- thetai*p.k*p.i*p.i.ij*(1-2*p.i+p.j.ij)+thetai*p.k*p.i.ij*p.j.ij*(p.i-p.j.ij+p.i.ij)
    iajd <- jdia <- -thetai*p.k*p.i.ij*(-1*p.j*p.j.ij+2*p.j*p.i-p.j.ij*p.j.ij+p.j.ij*p.i+p.i.ij*p.j.ij)
    #jakd <- kdja <- 2*thetaj*p.k*p.i*p.j
    jaid <- idja <- -thetaj*p.k*p.i.ij*(-1*p.j*p.j.ij+2*p.j*p.i-p.j.ij*p.j.ij+p.j.ij*p.i+p.i.ij*p.j.ij)
    jajd <- jdja <- -thetaj*p.k*p.i.ij*((p.j+p.j.ij)^2+p.j^2-p.j-p.i.ij*p.j.ij)
    #kdid <- idkd <- p.k*p.i*(1-2*p.i)
    #kdjd <- jdkd <- -2*p.k*p.i*p.j
    idjd <- jdid <- p.k*p.i.ij*(-1*p.j*p.j.ij+2*p.j*p.i-p.j.ij*p.j.ij+p.j.ij*p.i+p.i.ij*p.j.ij)
  }else if(identical(x,c(3,2,1))){
    kaka <- thetak^2*p.k*p.j*(1-2*p.k)
    jaja <- -thetaj^2*p.k*p.j*p.j.ji*(1-2*p.j+p.i.ji)-thetaj^2*p.k*p.j.ji*p.i.ji*(p.j-p.i.ji+p.j.ji)
    iaia <- thetai^2*p.k*p.j.ji*((p.i+p.i.ji)^2+p.i^2-p.i-p.j.ji*p.i.ji)
    #kdkd <- p.k*p.j*(1-2*p.k)
    jdjd <- -p.k*p.j*p.j.ji*(1-2*p.j+p.i.ji)-p.k*p.j.ji*p.i.ji*(p.j-p.i.ji+p.j.ji)
    idid <- p.k*p.j.ji*((p.i+p.i.ji)^2+p.i^2-p.i-p.j.ji*p.i.ji)
    kaja <- jaka <- thetak*thetaj*p.k*p.j*(1-2*p.j)
    kaia <- iaka <- -2*thetak*thetai*p.k*p.j*p.i
    #kakd <- kdka <- -thetak*p.k*p.j*(1-2*p.k)
    kajd <- jdka <- -thetak*p.k*p.j*(1-2*p.j)
    kaid <- idka <- 2*thetak*p.k*p.j*p.i
    jaia <- iaja <- thetaj*thetai*p.k*p.j.ji*(-1*p.i*p.i.ji+2*p.i*p.j-p.i.ji*p.i.ji+p.i.ji*p.j+p.j.ji*p.i.ji)
    #jakd <- kdja <- -thetaj*p.k*p.j*(1-2*p.j)
    jajd <- jdja <- thetaj*p.k*p.j*p.j.ji*(1-2*p.j+p.i.ji)+thetaj*p.k*p.j.ji*p.i.ji*(p.j-p.i.ji+p.j.ji)
    jaid <- idja <- -thetaj*p.k*p.j.ji*(-1*p.i*p.i.ji+2*p.i*p.j-p.i.ji*p.i.ji+p.i.ji*p.j+p.j.ji*p.i.ji)
    #iakd <- kdia <- 2*thetai*p.k*p.j*p.i
    iajd <- jdia <- -thetai*p.k*p.j.ji*(-1*p.i*p.i.ji+2*p.i*p.j-p.i.ji*p.i.ji+p.i.ji*p.j+p.j.ji*p.i.ji)
    iaid <- idia <- -thetai*p.k*p.j.ji*((p.i+p.i.ji)^2+p.i^2-p.i-p.j.ji*p.i.ji)
    #kdjd <- jdkd <- p.k*p.j*(1-2*p.j)
    #kdid <- idkd <- -2*p.k*p.j*p.i
    jdid <- idjd <- p.k*p.j.ji*(-1*p.i*p.i.ji+2*p.i*p.j-p.i.ji*p.i.ji+p.i.ji*p.j+p.j.ji*p.i.ji)
  }
  h <- diag(c(iaia,jaja,kaka,idid,jdjd))
  h[1,2] <- iaja
  h[1,3] <- iaka
  h[1,4] <- iaid
  h[1,5] <- iajd
  h[2,1] <- jaia
  h[2,3] <- jaka
  h[2,4] <- jaid
  h[2,5] <- jajd
  h[3,1] <- kaia
  h[3,2] <- kaja
  h[3,4] <- kaid
  h[3,5] <- kajd
  h[4,1] <- idia
  h[4,2] <- idja
  h[4,3] <- idka
  h[4,5] <- idjd
  h[5,1] <- jdia
  h[5,2] <- jdja
  h[5,3] <- jdka
  h[5,4] <- jdid
  h
}

# information matrix of block b

## 入参：以block为单位计算
## theta=作答者能力值，a=题目的a参数，d=（ a * b）


Icom_item <- function(theta,a,d){
  all.x <- list(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
  ret <- 0
  for(x in 1:length(all.x)){
    g <- grad_item(theta,a,d,x=all.x[[x]])
    ret <- ret + tcrossprod(g,g)/px(theta,a,d,xx=all.x[[x]])-hes_item(theta,a,d,x=all.x[[x]])
  }
  ret
}
Ess_item <- function(theta,a,d){
  all.x <- list(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
  ret <- 0
  for(x in 1:length(all.x)){
    g <- grad_item(theta,a,d,x=all.x[[x]])
    ret <- ret + tcrossprod(g,g)/px(theta,a,d,xx=all.x[[x]])
  }
  ret
}
SS_item <- function(theta,a,d,Y){
  all.x <- list(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
  ret <- 0
  
  g <- grad_item(theta,a,d,x=all.x[[Y]])
  ret <- tcrossprod(g,g)/px(theta,a,d,xx=all.x[[Y]])^2
  
  ret
}

# Ib(theta=c(0,0,0),a=c(1,2,3),d=c(0.5,1,2))


# test information matrix = sum of Ib
# this also returns SE of each dimension
# 给定theta计算相应的se，各block的信息量

## 入参：
## theta=作答者能力值，item.par=全卷题目参数表（见使用示例），BID=block-题目-维度对应关系矩阵


test.info_item <- function(theta,item.par,BID,Y){
  D <- max(BID$Block)*(max(BID$Item)*2-1) #number of dims
  complete.info.b <- diag(rep(0,D))
  for(b in unique(BID$Block)){
    items <- which(BID$Block==b)
    item <- c((b-1)*5+1,(b-1)*5+2,(b-1)*5+3,(b-1)*5+4,(b-1)*5+5)
    dim.b <- BID$Dim[items]  
    a.b <- item.par[items,1]
    d.b <- item.par[items,2]
    theta.b <- theta[dim.b]
    Y.b <- Y[b]
    complete.info.b[item,item] <- SS_item(theta.b,a.b,d.b,Y.b) 
  }
  complete.info.b
}

singleI <- function(theta,item.par,BID,Y){
  
  #sinI=matrix(0,max(BID$Block)*(max(BID$Item)*2-1),max(BID$Block)*(max(BID$Item)*2-1))
  sinI=foreach(i = seq_len(nrow(theta)),.inorder=FALSE,.combine='+',.export = c("test.info_item","Ib_item","px","grad_item","hes_item","Y"))%dopar%{
    test.info_item(theta[i,],item.par,BID,Y[i,])
  }
  sinI
}

SE_item_3 <- function(x,BID,Y){
  library(foreach)
  library(doParallel)
  nCPUcores = detectCores()
  if (nCPUcores < 3) {
    registerDoSEQ()
  }else{
    cl = makeCluster(1)
    registerDoParallel(cl)
  }
  
  all=foreach(i = seq_len(length(x$thetalist)),.inorder=FALSE,.combine='+',.export = c("foreach","%dopar%","singleI","test.info_item","Ib_item","px","hes_item","grad_item","Y"))%dopar%{
    singleI(x$thetalist[[i]],data.frame(c(x$a),c(x$d)),BID,Y)
  }
  meanI=all/length(x$thetalist)
  se <- matrix(sqrt(diag(solve(meanI))),nrow=max(BID$Block),ncol=5,byrow = TRUE)
  se.a <- se[,1:3]
  se.d <- se[,4:5]
  se.d <- cbind(se.d,sqrt(se.d[,1]^2+se.d[,2]^2))
  list(se.a=se.a,se.d=se.d)
}


