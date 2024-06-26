
px_2PL_2 <- function(theta,a,d,x=1){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)
  z <- exp(a*theta-d)          #1*3

  # p.j.ijk <- p.j.ikj <- p.j.jki <- p.j.jik <- p.j.kij <- p.j.kji
  p.j <- z[2]/sum(z)
  # p.i.ijk <- p.i.ikj <- p.i.jki <- p.i.jik <- p.i.kij <- p.i.kji
  p.i <- z[1]/sum(z)

  if(identical(x,1)){       #顺序为1,2,3时的一阶导，下方同理
    L=p.i
  }else if(identical(x,2)){
    L=p.j
  }
  L
}


singleI_2PL_2 <- function(theta,item.par,BID,Y,model,I,h=NULL){
  # information matrix of block b

  ## 入参：以block为单位计算
  ## theta=作答者能力值，a=题目的a参数，d=（ a * b）

  grad_FDM_item_2PL_2 <- function(theta,a,d,x=1,h){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)

    g.ia <- (px_2PL_2(theta,a+c(h,0),d,x)-px_2PL_2(theta,a,d,x))/h
    g.ja <- (px_2PL_2(theta,a+c(0,h),d,x)-px_2PL_2(theta,a,d,x))/h
    g.id <- (px_2PL_2(theta,a,d+c(h,-h),x)-px_2PL_2(theta,a,d,x))/h

    c(g.ia,g.ja,g.id)
  }
  grad_CDM_item_2PL_2 <- function(theta,a,d,x=1,h){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)

    g.ia <- (px_2PL_2(theta,a+c(h,0),d,x)-px_2PL_2(theta,a-c(h,0),d,x))/2*h
    g.ja <- (px_2PL_2(theta,a+c(0,h),d,x)-px_2PL_2(theta,a-c(0,h),d,x))/2*h
    g.id <- (px_2PL_2(theta,a,d+c(h,-h),x)-px_2PL_2(theta,a,d-c(h,-h),x))/2*h

    c(g.ia,g.ja,g.id)
  }
  grad_RES_item_2PL_2 <- function(theta,a,d,x=1,h){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)

    g.ia <- (px_2PL_2(theta,a-c(2*h,0),d,x)+8*px_2PL_2(theta,a+c(h,0),d,x)
             -px_2PL_2(theta,a+c(2*h,0),d,x)-8*px_2PL_2(theta,a-c(h,0),d,x))/12*h
    g.ja <- (px_2PL_2(theta,a-c(0,2*h),d,x)+8*px_2PL_2(theta,a+c(0,h),d,x)
             -px_2PL_2(theta,a+c(0,2*h),d,x)-8*px_2PL_2(theta,a-c(0,h),d,x))/12*h
    g.id <- (px_2PL_2(theta,a,d-c(2*h,-2*h),x)+8*px_2PL_2(theta,a,d+c(h,-h),x)
             -px_2PL_2(theta,a,d+c(2*h,-2*h),x)-8*px_2PL_2(theta,a,d-c(h,-h),x))/12*h
    c(g.ia,g.ja,g.id)
  }


  grad_item_2PL_2 <- function(theta,a,d,x=1){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)
    z <- exp(a*theta-d)          #1*3
    thetai <- theta[1]
    thetaj <- theta[2]

    # p.j.ijk <- p.j.ikj <- p.j.jki <- p.j.jik <- p.j.kij <- p.j.kji
    p.j <- z[2]/sum(z)
    # p.i.ijk <- p.i.ikj <- p.i.jki <- p.i.jik <- p.i.kij <- p.i.kji
    p.i <- z[1]/sum(z)

    if(identical(x,1)){
      g.ia <- thetai*p.i*p.j
      g.ja <- -thetaj*p.i*p.j
      g.id <- -2*p.i*p.j
    }else if(identical(x,2)){
      g.ja <- thetaj*p.i*p.j
      g.ia <- -thetai*p.j*p.i
      g.id <- 2*p.i*p.j
    }

    c(g.ia,g.ja,g.id)
  }



  FDMEss_item_2PL_2 <- function(theta,a,d,Y,h=1e-5){
    all.x <- list(1,2)
    g <- grad_FDM_item_2PL_2(theta,a,d,x=all.x[[Y]],h=h)
    ret <- tcrossprod(g,g)/px_2PL_2(theta,a,d,xx=all.x[[Y]])^2

    ret
  }
  RESEss_item_2PL_2 <- function(theta,a,d,Y,h=1){
    all.x <- list(1,2)
    g <- grad_RES_item_2PL_2(theta,a,d,x=all.x[[Y]],h=h)
    ret <- tcrossprod(g,g)/px_2PL_2(theta,a,d,xx=all.x[[Y]])^2

    ret
  }
  CDMEss_item_2PL_2 <- function(theta,a,d,Y,h=1){
    all.x <- list(1,2)
    g <- grad_CDM_item_2PL_2(theta,a,d,x=all.x[[Y]],h=h)
    ret <- tcrossprod(g,g)/px_2PL_2(theta,a,d,xx=all.x[[Y]])^2

    ret
  }
  XPDEss_item_2PL_2 <- function(theta,a,d,Y){
    all.x <- list(1,2)
    g <- grad_item_2PL_2(theta,a,d,x=all.x[[Y]])
    ret <- tcrossprod(g,g)/px_2PL_2(theta,a,d,xx=all.x[[Y]])^2

    ret
  }


  FDMSS_item_2PL_2 <- function(theta,a,d,Y,h=1e-5){
    all.x <- list(1,2)
    g <- grad_FDM_item_2PL_2(theta,a,d,x=all.x[[Y]],h=h)
    g
  }
  RESSS_item_2PL_2 <- function(theta,a,d,Y,h=1){
    all.x <- list(1,2)
    g <- grad_RES_item_2PL_2(theta,a,d,x=all.x[[Y]],h=h)
    g
  }
  CDMSS_item_2PL_2 <- function(theta,a,d,Y,h=1){
    all.x <- list(1,2)
    g <- grad_CDM_item_2PL_2(theta,a,d,x=all.x[[Y]],h=h)
    g
  }
  XPDSS_item_2PL_2 <- function(theta,a,d,Y){
    all.x <- list(1,2)
    g <- grad_item_2PL_2(theta,a,d,x=all.x[[Y]])
    g
  }


  # Ib(theta=c(0,0,0),a=c(1,2,3),d=c(0.5,1,2))


  # test information matrix = sum of Ib
  # this also returns SE of each dimension
  # 给定theta计算相应的se，各block的信息量

  ## 入参：
  ## theta=作答者能力值，item.par=全卷题目参数表（见使用示例），BID=block-题目-维度对应关系矩阵



  itemI_2PL_2 <- function(theta,item.par,BID,Y,I,h=NULL){
    D <- max(BID$Block)*(max(BID$Item)*2-1) #number of dims
    itemIEss <- diag(rep(0,D))
    itemISS <- rep(0,D)

    for(b in unique(BID$Block)){
      items <- which(BID$Block==b)
      item <- c((b-1)*3+1,(b-1)*3+2,(b-1)*3+3)
      dim.b <- BID$Dim[items]
      a.b <- item.par[items,1]
      d.b <- item.par[items,2]
      theta.b <- theta[dim.b]
      Y.b <- Y[b]
      if(I=="XPD"){
        itemIEss[item,item] <- XPDEss_item_2PL_2(theta.b,a.b,d.b,Y.b)
        itemISS[item] <- XPDSS_item_2PL_2(theta.b,a.b,d.b,Y.b)
      }else if(I=="FDM"){
        if(is.null(h)){
          itemIEss[item,item] <- FDMEss_item_2PL_2(theta.b,a.b,d.b,Y.b)
          itemISS[item] <- FDMSS_item_2PL_2(theta.b,a.b,d.b,Y.b)
        }else{
          itemIEss[item,item] <- FDMEss_item_2PL_2(theta.b,a.b,d.b,Y.b,h=h)
          itemISS[item] <- FDMSS_item_2PL_2(theta.b,a.b,d.b,Y.b,h=h)
        }
      }else if(I=="CDM"){
        if(is.null(h)){
          itemIEss[item,item] <- CDMEss_item_2PL_2(theta.b,a.b,d.b,Y.b)
          itemISS[item] <- CDMSS_item_2PL_2(theta.b,a.b,d.b,Y.b)
        }else{
          itemIEss[item,item] <- CDMEss_item_2PL_2(theta.b,a.b,d.b,Y.b,h=h)
          itemISS[item] <- CDMSS_item_2PL_2(theta.b,a.b,d.b,Y.b,h=h)
        }
      }else if(I=="RES"){
        if(is.null(h)){
          itemIEss[item,item] <- RESEss_item_2PL_2(theta.b,a.b,d.b,Y.b)
          itemISS[item] <- RESSS_item_2PL_2(theta.b,a.b,d.b,Y.b)
        }else{
          itemIEss[item,item] <- RESEss_item_2PL_2(theta.b,a.b,d.b,Y.b,h=h)
          itemISS[item] <- RESSS_item_2PL_2(theta.b,a.b,d.b,Y.b,h=h)
        }
      }
    }
    list(itemIEss=itemIEss,itemISS=itemISS)
  }
  sinI_matrix=foreach(i = seq_len(nrow(theta)),.inorder=FALSE,.combine='+',.export = c("foreach","%dopar%","px_2PL_2"))%dopar%{
                                                                                         itemI_2PL_2(theta[i,],item.par,BID,Y[i,],I,h=h)$itemIEss
                                                                                       }
  sinI_vector=foreach(i = seq_len(nrow(theta)),.inorder=FALSE,.combine='+',.export = c("foreach","%dopar%","px_2PL_2"))%dopar%{
                                                                                         itemI_2PL_2(theta[i,],item.par,BID,Y[i,],I,h=h)$itemISS
                                                                                       }
  list(sinI_matrix=sinI_matrix,sinI_vector=sinI_vector)
}

itemSE_2PL_2 <- function(x,BID,Y,model,SE,IN,theta,h=NULL){


  matrix=singleI_2PL_2(theta=theta,item.par=data.frame(c(x$a),c(x$d)),BID=BID,Y=Y,I=SE,model=model,h=h)$sinI_matrix
  vector=singleI_2PL_2(theta=theta,item.par=data.frame(c(x$a),c(x$d)),BID=BID,Y=Y,I=SE,model=model,h=h)$sinI_vector
  I=matrix-tcrossprod(vector,vector)/nrow(theta)

  se <- matrix(sqrt(diag(solve(I))),nrow=max(BID$Block),ncol=3,byrow = TRUE)
  se.a <- se[,1:2]
  se.d <- se[,3]
  se.d <- cbind(se.d,se.d)
  list(SE.a=se.a,SE.d=se.d,SE.sigm=NULL)
}


