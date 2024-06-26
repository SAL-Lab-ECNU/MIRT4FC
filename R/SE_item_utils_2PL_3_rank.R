

singleI_2PL_3_rank <- function(theta,item.par,BID,Y,model,I,h=NULL){
  grad_FDM_item_2PL_3_rank <- function(theta,a,d,x=c(1,2,3),h){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)

    g.ia <- (px_2PL_3_rank(theta,a+c(h,0,0),d,x)-px_2PL_3_rank(theta,a,d,x))/h
    g.ja <- (px_2PL_3_rank(theta,a+c(0,h,0),d,x)-px_2PL_3_rank(theta,a,d,x))/h
    g.ka <- (px_2PL_3_rank(theta,a+c(0,0,h),d,x)-px_2PL_3_rank(theta,a,d,x))/h
    g.id <- (px_2PL_3_rank(theta,a,d+c(h,0,-h),x)-px_2PL_3_rank(theta,a,d,x))/h
    g.jd <- (px_2PL_3_rank(theta,a,d+c(0,h,-h),x)-px_2PL_3_rank(theta,a,d,x))/h
    c(g.ia,g.ja,g.ka,g.id,g.jd)
  }
  grad_CDM_item_2PL_3_rank <- function(theta,a,d,x=c(1,2,3),h){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)

    g.ia <- (px_2PL_3_rank(theta,a+c(h,0,0),d,x)-px_2PL_3_rank(theta,a-c(h,0,0),d,x))/2*h
    g.ja <- (px_2PL_3_rank(theta,a+c(0,h,0),d,x)-px_2PL_3_rank(theta,a-c(0,h,0),d,x))/2*h
    g.ka <- (px_2PL_3_rank(theta,a+c(0,0,h),d,x)-px_2PL_3_rank(theta,a-c(0,0,h),d,x))/2*h
    g.id <- (px_2PL_3_rank(theta,a,d+c(h,0,-h),x)-px_2PL_3_rank(theta,a,d-c(h,0,-h),x))/2*h
    g.jd <- (px_2PL_3_rank(theta,a,d+c(0,h,-h),x)-px_2PL_3_rank(theta,a,d-c(0,h,-h),x))/2*h

    c(g.ia,g.ja,g.ka,g.id,g.jd)
  }
  grad_RES_item_2PL_3_rank <- function(theta,a,d,x=c(1,2,3),h){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)

    g.ia <- (px_2PL_3_rank(theta,a-c(2*h,0,0),d,x)+8*px_2PL_3_rank(theta,a+c(h,0,0),d,x)
             -px_2PL_3_rank(theta,a+c(2*h,0,0),d,x)-8*px_2PL_3_rank(theta,a-c(h,0,0),d,x))/12*h
    g.ja <- (px_2PL_3_rank(theta,a-c(0,2*h,0),d,x)+8*px_2PL_3_rank(theta,a+c(0,h,0),d,x)
             -px_2PL_3_rank(theta,a+c(0,2*h,0),d,x)-8*px_2PL_3_rank(theta,a-c(0,h,0),d,x))/12*h
    g.ka <- (px_2PL_3_rank(theta,a-c(0,0,2*h),d,x)+8*px_2PL_3_rank(theta,a+c(0,0,h),d,x)
             -px_2PL_3_rank(theta,a+c(0,0,2*h),d,x)-8*px_2PL_3_rank(theta,a-c(0,0,h),d,x))/12*h
    g.id <- (px_2PL_3_rank(theta,a,d-c(2*h,0,-2*h),x)+8*px_2PL_3_rank(theta,a,d+c(h,0,-h),x)
             -px_2PL_3_rank(theta,a,d+c(2*h,0,-2*h),x)-8*px_2PL_3_rank(theta,a,d-c(h,0,-h),x))/12*h
    g.jd <- (px_2PL_3_rank(theta,a,d-c(0,2*h,-2*h),x)+8*px_2PL_3_rank(theta,a,d+c(0,h,-h),x)
             -px_2PL_3_rank(theta,a,d+c(0,2*h,-2*h),x)-8*px_2PL_3_rank(theta,a,d-c(0,h,-h),x))/12*h

    c(g.ia,g.ja,g.ka,g.id,g.jd)
  }


  grad_item_2PL_3_rank <- function(theta,a,d,x=c(1,2,3)){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)
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
      g.id <- -p.i*p.j.jk*(p.j+2*p.k+p.k.jk)
      g.jd <- p.i*p.j.jk*(p.j-p.k-2*p.k.jk)
      #g.kd <- p.i*p.j.jk*(p.k+p.k.jk)
    }else if(identical(x,c(2,3,1))){
      g.ja <- thetaj*p.j*p.k
      g.ka <- thetak*p.j*p.k.ki*(p.i.ki-p.k)
      g.ia <- -thetai*p.j*p.k.ki*(p.i+p.i.ki)
      g.jd <- -p.j*p.k.ik*(p.i+2*p.k-p.i.ik)
      #g.kd <- -p.j*p.k.ki*(p.i.ki-p.k)
      g.id <- p.j*p.k.ki*(p.i-p.k+2*p.i.ki)
    }else if(identical(x,c(1,3,2))){
      g.ia <- thetai*p.i*p.k
      g.ka <- thetak*p.i*p.k.kj*(p.j.kj-p.k)
      g.ja <- -thetaj*p.i*p.k.kj*(p.j+p.j.kj)
      g.id <- -p.i*p.k.jk*(p.j+2*p.k-p.j.jk)
      #g.kd <- -p.i*p.k.kj*(p.j.kj-p.k)
      g.jd <- p.i*p.k.kj*(p.j-p.k+2*p.j.kj)
    }else if(identical(x,c(2,1,3))){
      g.ja <- thetaj*p.j*p.i
      g.ia <- thetai*p.j*p.i.ik*(p.k.ik-p.i)
      g.ka <- -thetak*p.j*p.i.ik*(p.k+p.k.ik)
      g.jd <- -p.j*p.i.ik*(p.i+2*p.k+p.k.ik)
      g.id <- p.j*p.i.ik*(p.i-p.k-2*p.k.ik)
      #g.kd <- p.j*p.i.ik*(p.k+p.k.ik)
    }else if(identical(x,c(3,1,2))){
      g.ka <- thetak*p.k*p.i
      g.ia <- thetai*p.k*p.i.ij*(p.j.ij-p.i)
      g.ja <- -thetaj*p.k*p.i.ij*(p.j+p.j.ij)
      #g.kd <- -p.k*p.i
      g.id <- p.k*p.i.ij*(2*p.i+p.j-p.j.ij)
      g.jd <- p.k*p.i.ij*(2*p.j+p.i+p.j.ij)
    }else if(identical(x,c(3,2,1))){
      g.ka <- thetak*p.k*p.j
      g.ja <- thetaj*p.k*p.j.ji*(p.i.ji-p.j)
      g.ia <- -thetai*p.k*p.j.ji*(p.i+p.i.ji)
      #g.kd <- -p.k*p.j
      g.jd <- p.k*p.j.ij*(2*p.j+p.i-p.i.ij)
      g.id <- p.k*p.j.ij*(2*p.i+p.j+p.i.ij)
    }

    c(g.ia,g.ja,g.ka,g.id,g.jd)
  }

  XPDEss_item_2PL_3_rank <- function(theta,a,d,Y){
    all.x <- list(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
    g <- grad_item_2PL_3_rank(theta,a,d,x=all.x[[Y]])
    ret <- tcrossprod(g,g)/px_2PL_3_rank(theta,a,d,xx=all.x[[Y]])^2

    ret
  }
  FDMEss_item_2PL_3_rank <- function(theta,a,d,Y,h=1e-5){
    all.x <- list(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
    g <- grad_FDM_item_2PL_3_rank(theta,a,d,x=all.x[[Y]],h=h)
    ret <- tcrossprod(g,g)/px_2PL_3_rank(theta,a,d,xx=all.x[[Y]])^2

    ret
  }
  RESEss_item_2PL_3_rank <- function(theta,a,d,Y,h=1){
    all.x <- list(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
    g <- grad_RES_item_2PL_3_rank(theta,a,d,x=all.x[[Y]],h=h)
    ret <- tcrossprod(g,g)/px_2PL_3_rank(theta,a,d,xx=all.x[[Y]])^2

    ret
  }
  CDMEss_item_2PL_3_rank <- function(theta,a,d,Y,h=1){
    all.x <- list(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
    g <- grad_CDM_item_2PL_3_rank(theta,a,d,x=all.x[[Y]],h=h)
    ret <- tcrossprod(g,g)/px_2PL_3_rank(theta,a,d,xx=all.x[[Y]])^2

    ret
  }



  FDMSS_item_2PL_3_rank <- function(theta,a,d,Y,h=1e-5){
    all.x <- list(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
    g <- grad_FDM_item_2PL_3_rank(theta,a,d,x=all.x[[Y]],h=h)
    g
  }
  RESSS_item_2PL_3_rank <- function(theta,a,d,Y,h=1){
    all.x <- list(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
    g <- grad_RES_item_2PL_3_rank(theta,a,d,x=all.x[[Y]],h=h)
    g
  }
  CDMSS_item_2PL_3_rank <- function(theta,a,d,Y,h=1){
    all.x <- list(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
    g <- grad_CDM_item_2PL_3_rank(theta,a,d,x=all.x[[Y]],h=h)
    g
  }
  XPDSS_item_2PL_3_rank <- function(theta,a,d,Y){
    all.x <- list(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
    g <- grad_item_2PL_3_rank(theta,a,d,x=all.x[[Y]])
    g
  }


  itemI_2PL_3_rank <- function(theta,item.par,BID,Y,I,h=NULL){
    D <- max(BID$Block)*(max(BID$Item)*2-1) #number of dims
    itemIEss <- diag(rep(0,D))
    itemISS <- rep(0,D)
    for(b in unique(BID$Block)){
      items <- which(BID$Block==b)
      item <- c((b-1)*5+1,(b-1)*5+2,(b-1)*5+3,(b-1)*5+4,(b-1)*5+5)
      dim.b <- BID$Dim[items]
      a.b <- item.par[items,1]
      d.b <- item.par[items,2]
      theta.b <- theta[dim.b]
      Y.b <- Y[b]
      if(I=="XPD"){
        itemIEss[item,item] <- XPDEss_item_2PL_3_rank(theta.b,a.b,d.b,Y.b)
        itemISS[item] <- XPDSS_item_2PL_3_rank(theta.b,a.b,d.b,Y.b)
      }else if(I=="FDM"){
        if(is.null(h)){
          itemIEss[item,item] <- FDMEss_item_2PL_3_rank(theta.b,a.b,d.b,Y.b)
          itemISS[item] <- FDMSS_item_2PL_3_rank(theta.b,a.b,d.b,Y.b)
        }else{
          itemIEss[item,item] <- FDMEss_item_2PL_3_rank(theta.b,a.b,d.b,Y.b,h=h)
          itemISS[item] <- FDMSS_item_2PL_3_rank(theta.b,a.b,d.b,Y.b,h=h)
        }
      }else if(I=="CDM"){
        if(is.null(h)){
          itemIEss[item,item] <- CDMEss_item_2PL_3_rank(theta.b,a.b,d.b,Y.b)
          itemISS[item] <- CDMSS_item_2PL_3_rank(theta.b,a.b,d.b,Y.b)
        }else{
          itemIEss[item,item] <- CDMEss_item_2PL_3_rank(theta.b,a.b,d.b,Y.b,h=h)
          itemISS[item] <- CDMSS_item_2PL_3_rank(theta.b,a.b,d.b,Y.b,h=h)
        }
      }else if(I=="RES"){
        if(is.null(h)){
          itemIEss[item,item] <- RESEss_item_2PL_3_rank(theta.b,a.b,d.b,Y.b)
          itemISS[item] <- RESSS_item_2PL_3_rank(theta.b,a.b,d.b,Y.b)
        }else{
          itemIEss[item,item] <- RESEss_item_2PL_3_rank(theta.b,a.b,d.b,Y.b,h=h)
          itemISS[item] <- RESSS_item_2PL_3_rank(theta.b,a.b,d.b,Y.b,h=h)
        }
      }
    }
    list(itemIEss=itemIEss,itemISS=itemISS)
  }


  # information matrix of block b

  ## 入参：以block为单位计算
  ## theta=作答者能力值，a=题目的a参数，d=（ a * b）





  sinI_matrix=foreach(i = seq_len(nrow(theta)),.inorder=FALSE,.combine='+',.export = c("foreach","%dopar%","px_2PL_3_rank"))%dopar%{
    itemI_2PL_3_rank(theta[i,],item.par,BID,Y[i,],I,h=h)$itemIEss
  }
  sinI_vector=foreach(i = seq_len(nrow(theta)),.inorder=FALSE,.combine='+',.export = c("foreach","%dopar%","px_2PL_3_rank"))%dopar%{
    itemI_2PL_3_rank(theta[i,],item.par,BID,Y[i,],I,h=h)$itemISS
  }
  list(sinI_matrix=sinI_matrix,sinI_vector=sinI_vector)
}

itemSE_2PL_3_rank <- function(x,BID,Y,model,SE,IN,theta,h=NULL){


  matrix=singleI_2PL_3_rank(theta=theta,item.par=data.frame(c(x$a),c(x$d)),BID=BID,Y=Y,I=SE,model=model,h=h)$sinI_matrix
  vector=singleI_2PL_3_rank(theta=theta,item.par=data.frame(c(x$a),c(x$d)),BID=BID,Y=Y,I=SE,model=model,h=h)$sinI_vector
  I=matrix-tcrossprod(vector,vector)/nrow(theta)

  se <- matrix(sqrt(diag(solve(I))),nrow=max(BID$Block),ncol=5,byrow = TRUE)
  se.a <- se[,1:3]
  se.d <- se[,4:5]
  se.d <- cbind(se.d,sqrt(se.d[,1]^2+se.d[,2]^2))
  list(SE.a=se.a,SE.d=se.d,SE.sigm=NULL)
}

