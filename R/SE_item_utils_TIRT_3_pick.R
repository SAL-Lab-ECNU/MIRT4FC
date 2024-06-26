

singleI_TIRT_3_pick <- function(theta,item.par,BID,Y,model,I,h=NULL){
  grad_FDM_item_TIRT_3_pick <- function(theta,a,d,x=1,h){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)

    g.ia <- (px_TIRT_3_pick(theta,a+c(h,0,0),d,x)-px_TIRT_3_pick(theta,a,d,x))/h
    g.ja <- (px_TIRT_3_pick(theta,a+c(0,h,0),d,x)-px_TIRT_3_pick(theta,a,d,x))/h
    g.ka <- (px_TIRT_3_pick(theta,a+c(0,0,h),d,x)-px_TIRT_3_pick(theta,a,d,x))/h
    g.id <- (px_TIRT_3_pick(theta,a,d+c(h,0,-h),x)-px_TIRT_3_pick(theta,a,d,x))/h
    g.jd <- (px_TIRT_3_pick(theta,a,d+c(0,h,-h),x)-px_TIRT_3_pick(theta,a,d,x))/h
    #g.kd <- (px_TIRT_3_pick(theta,a,d+c(0,0,h),x)-px_TIRT_3_pick(theta,a,d,x))/h
    c(g.ia,g.ja,g.ka,g.id,g.jd)
  }
  grad_CDM_item_TIRT_3_pick <- function(theta,a,d,x=1,h){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)

    g.ia <- (px_TIRT_3_pick(theta,a+c(h,0,0),d,x)-px_TIRT_3_pick(theta,a-c(h,0,0),d,x))/2*h
    g.ja <- (px_TIRT_3_pick(theta,a+c(0,h,0),d,x)-px_TIRT_3_pick(theta,a-c(0,h,0),d,x))/2*h
    g.ka <- (px_TIRT_3_pick(theta,a+c(0,0,h),d,x)-px_TIRT_3_pick(theta,a-c(0,0,h),d,x))/2*h
    g.id <- (px_TIRT_3_pick(theta,a,d+c(h,0,-h),x)-px_TIRT_3_pick(theta,a,d-c(h,0,-h),x))/2*h
    g.jd <- (px_TIRT_3_pick(theta,a,d+c(0,h,-h),x)-px_TIRT_3_pick(theta,a,d-c(0,h,-h),x))/2*h
    #g.kd <- (px_TIRT_3_pick(theta,a,d+c(0,0,h),x)-px_TIRT_3_pick(theta,a,d-c(0,0,h),x))/2*h
    c(g.ia,g.ja,g.ka,g.id,g.jd)
  }
  grad_RES_item_TIRT_3_pick <- function(theta,a,d,x=1,h){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)

    g.ia <- (px_TIRT_3_pick(theta,a-c(2*h,0,0),d,x)+8*px_TIRT_3_pick(theta,a+c(h,0,0),d,x)
             -px_TIRT_3_pick(theta,a+c(2*h,0,0),d,x)-8*px_TIRT_3_pick(theta,a-c(h,0,0),d,x))/12*h
    g.ja <- (px_TIRT_3_pick(theta,a-c(0,2*h,0),d,x)+8*px_TIRT_3_pick(theta,a+c(0,h,0),d,x)
             -px_TIRT_3_pick(theta,a+c(0,2*h,0),d,x)-8*px_TIRT_3_pick(theta,a-c(0,h,0),d,x))/12*h
    g.ka <- (px_TIRT_3_pick(theta,a-c(0,0,2*h),d,x)+8*px_TIRT_3_pick(theta,a+c(0,0,h),d,x)
             -px_TIRT_3_pick(theta,a+c(0,0,2*h),d,x)-8*px_TIRT_3_pick(theta,a-c(0,0,h),d,x))/12*h
    g.id <- (px_TIRT_3_pick(theta,a,d-c(2*h,0,-2*h),x)+8*px_TIRT_3_pick(theta,a,d+c(h,0,-h),x)
             -px_TIRT_3_pick(theta,a,d+c(2*h,0,-2*h),x)-8*px_TIRT_3_pick(theta,a,d-c(h,0,-h),x))/12*h
    g.jd <- (px_TIRT_3_pick(theta,a,d-c(0,2*h,-2*h),x)+8*px_TIRT_3_pick(theta,a,d+c(0,h,-h),x)
             -px_TIRT_3_pick(theta,a,d+c(0,2*h,-2*h),x)-8*px_TIRT_3_pick(theta,a,d-c(0,h,-h),x))/12*h
    #g.kd <- (px_TIRT_3_pick(theta,a,d-c(0,0,2*h),x)+8*px_TIRT_3_pick(theta,a,d+c(0,0,h),x)
    #         -px_TIRT_3_pick(theta,a,d+c(0,0,2*h),x)-8*px_TIRT_3_pick(theta,a,d-c(0,0,h),x))/12*h
    c(g.ia,g.ja,g.ka,g.id,g.jd)
  }

  grad_item_TIRT_3_pick <- function(theta,a,d,x=1){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)
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
    thetai <- theta[1]
    thetaj <- theta[2]
    thetak <- theta[3]

    if(identical(x,1)){       #顺序为1,2,3时的一阶导，下方同理
      g.ia <- (thetai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13+p12*thetai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))/pall-
        pa/pall^2*((thetai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13+p12*thetai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+(-thetai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23)+(-thetai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32))
      g.ja <- (-thetaj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13)/pall-
        pa/pall^2*((-thetaj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13)+(thetaj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p23+p21*thetaj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+(-thetaj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)*p31))
      g.ka <- (-p12*thetak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))/pall-
        pa/pall^2*((-p12*thetak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+(-p21*thetak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+(p32*thetak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)+p31*thetak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)))
      g.id <- (-1/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13-2*p12/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))/pall-
        pa/pall^2*((-1/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13-2*p12/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+(1/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23-p21/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+(2/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32+p31/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)))
      g.jd <- (1/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13-p12/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))/pall-
        pa/pall^2*((1/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13-p12/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+(-1/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23-2*p21/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+(1/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32+2*p31/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)))

    }else if(identical(x,2)){
      g.ia <- (-thetai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23)/pall-
        pb/pall^2*((thetai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13+p12*thetai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+(-thetai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23)+(-thetai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32))
      g.ja <- (thetaj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p23+p21*thetaj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))/pall-
        pb/pall^2*((-thetaj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13)+(thetaj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p23+p21*thetaj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+(-thetaj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)*p31))
      g.ka <- (-p21*thetak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))/pall-
        pb/pall^2*((-p12*thetak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+(-p21*thetak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+(p32*thetak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)+p31*thetak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)))
      g.id <- (1/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23-p21/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))/pall-
        pb/pall^2*((-1/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13-2*p12/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+(1/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23-p21/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+(2/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32+p31/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)))
      g.jd <- (-1/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23-2*p21/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))/pall-
        pb/pall^2*((1/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13-p12/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+(-1/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23-2*p21/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+(1/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32+2*p31/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)))

    }else if(identical(x,3)){
      g.ia <- (-thetai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32)/pall-
        pc/pall^2*((thetai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13+p12*thetai/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+(-thetai/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23)+(-thetai/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32))
      g.ja <- (-thetaj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)*p31)/pall-
        pc/pall^2*((-thetaj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13)+(thetaj/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p23+p21*thetaj/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+(-thetaj/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)*p31))
      g.ka <- (p32*thetak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)+p31*thetak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))/pall-
        pc/pall^2*((-p12*thetak/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+(-p21*thetak/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+(p32*thetak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)+p31*thetak/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)))
      g.id <- (2/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32+p31/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))/pall-
        pc/pall^2*((-1/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13-2*p12/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+(1/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23-p21/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+(2/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32+p31/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)))
      g.jd <- (1/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32+2*p31/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2))/pall-
        pc/pall^2*((1/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[2]-d[1]+d[2])^2/2)*p13-p12/sqrt(2*pi)*exp(-(alphatheta[1]-alphatheta[3]-d[1]+d[3])^2/2))+(-1/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[1]-d[2]+d[1])^2/2)*p23-2*p21/sqrt(2*pi)*exp(-(alphatheta[2]-alphatheta[3]-d[2]+d[3])^2/2))+(1/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[1]-d[3]+d[1])^2/2)*p32+2*p31/sqrt(2*pi)*exp(-(alphatheta[3]-alphatheta[2]-d[3]+d[2])^2/2)))
    }

    c(g.ia,g.ja,g.ka,g.id,g.jd)
  }


  XPDEss_item_TIRT_3_pick <- function(theta,a,d,Y){
    all.x <- list(1,2,3)
    g <- grad_item_TIRT_3_pick(theta,a,d,x=all.x[[Y]])
    ret <- tcrossprod(g,g)/px_TIRT_3_pick(theta,a,d,xx=all.x[[Y]])^2

    ret
  }
  FDMEss_item_TIRT_3_pick <- function(theta,a,d,Y,h=1e-5){
    all.x <- list(1,2,3)
    g <- grad_FDM_item_TIRT_3_pick(theta,a,d,x=all.x[[Y]],h=h)
    ret <- tcrossprod(g,g)/px_TIRT_3_pick(theta,a,d,xx=all.x[[Y]])^2

    ret
  }
  RESEss_item_TIRT_3_pick <- function(theta,a,d,Y,h=1){
    all.x <- list(1,2,3)
    g <- grad_RES_item_TIRT_3_pick(theta,a,d,x=all.x[[Y]],h=h)
    ret <- tcrossprod(g,g)/px_TIRT_3_pick(theta,a,d,xx=all.x[[Y]])^2

    ret
  }
  CDMEss_item_TIRT_3_pick <- function(theta,a,d,Y,h=1){
    all.x <- list(1,2,3)
    g <- grad_CDM_item_TIRT_3_pick(theta,a,d,x=all.x[[Y]],h=h)
    ret <- tcrossprod(g,g)/px_TIRT_3_pick(theta,a,d,xx=all.x[[Y]])^2

    ret
  }



  FDMSS_item_TIRT_3_pick <- function(theta,a,d,Y,h=1e-5){
    all.x <- list(1,2,3)
    g <- grad_FDM_item_TIRT_3_pick(theta,a,d,x=all.x[[Y]],h=h)
    g
  }
  RESSS_item_TIRT_3_pick <- function(theta,a,d,Y,h=1){
    all.x <- list(1,2,3)
    g <- grad_RES_item_TIRT_3_pick(theta,a,d,x=all.x[[Y]],h=h)
    g
  }
  CDMSS_item_TIRT_3_pick <- function(theta,a,d,Y,h=1){
    all.x <- list(1,2,3)
    g <- grad_CDM_item_TIRT_3_pick(theta,a,d,x=all.x[[Y]],h=h)
    g
  }
  XPDSS_item_TIRT_3_pick <- function(theta,a,d,Y){
    all.x <- list(1,2,3)
    g <- grad_item_TIRT_3_pick(theta,a,d,x=all.x[[Y]])
    g
  }


  itemI_TIRT_3_pick <- function(theta,item.par,BID,Y,I,h=NULL){
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
        itemIEss[item,item] <- XPDEss_item_TIRT_3_pick(theta.b,a.b,d.b,Y.b)
        itemISS[item] <- XPDSS_item_TIRT_3_pick(theta.b,a.b,d.b,Y.b)
      }else if(I=="FDM"){
        if(is.null(h)){
          itemIEss[item,item] <- FDMEss_item_TIRT_3_pick(theta.b,a.b,d.b,Y.b)
          itemISS[item] <- FDMSS_item_TIRT_3_pick(theta.b,a.b,d.b,Y.b)
        }else{
          itemIEss[item,item] <- FDMEss_item_TIRT_3_pick(theta.b,a.b,d.b,Y.b,h=h)
          itemISS[item] <- FDMSS_item_TIRT_3_pick(theta.b,a.b,d.b,Y.b,h=h)
        }
      }else if(I=="CDM"){
        if(is.null(h)){
          itemIEss[item,item] <- CDMEss_item_TIRT_3_pick(theta.b,a.b,d.b,Y.b)
          itemISS[item] <- CDMSS_item_TIRT_3_pick(theta.b,a.b,d.b,Y.b)
        }else{
          itemIEss[item,item] <- CDMEss_item_TIRT_3_pick(theta.b,a.b,d.b,Y.b,h=h)
          itemISS[item] <- CDMSS_item_TIRT_3_pick(theta.b,a.b,d.b,Y.b,h=h)
        }
      }else if(I=="RES"){
        if(is.null(h)){
          itemIEss[item,item] <- RESEss_item_TIRT_3_pick(theta.b,a.b,d.b,Y.b)
          itemISS[item] <- RESSS_item_TIRT_3_pick(theta.b,a.b,d.b,Y.b)
        }else{
          itemIEss[item,item] <- RESEss_item_TIRT_3_pick(theta.b,a.b,d.b,Y.b,h=h)
          itemISS[item] <- RESSS_item_TIRT_3_pick(theta.b,a.b,d.b,Y.b,h=h)
        }
      }
    }
    list(itemIEss=itemIEss,itemISS=itemISS)
  }


  # information matrix of block b

  ## 入参：以block为单位计算
  ## theta=作答者能力值，a=题目的a参数，d=（ a * b）





  sinI_matrix=foreach(i = seq_len(nrow(theta)),.inorder=FALSE,.combine='+',.export = c("px_TIRT_3_pick","foreach","%dopar%"))%dopar%{
                                                                                         itemI_TIRT_3_pick(theta[i,],item.par,BID,Y[i,],I,h=h)$itemIEss
                                                                                       }
  sinI_vector=foreach(i = seq_len(nrow(theta)),.inorder=FALSE,.combine='+',.export = c("px_TIRT_3_pick","foreach","%dopar%"))%dopar%{
                                                                                         itemI_TIRT_3_pick(theta[i,],item.par,BID,Y[i,],I,h=h)$itemISS
                                                                                       }
  list(sinI_matrix=sinI_matrix,sinI_vector=sinI_vector)
}

itemSE_TIRT_3_pick <- function(x,BID,Y,model,SE,IN,theta,h=NULL){


  matrix=singleI_TIRT_3_pick(theta=theta,item.par=data.frame(c(x$a),c(x$d)),BID=BID,Y=Y,I=SE,model=model,h=h)$sinI_matrix
  vector=singleI_TIRT_3_pick(theta=theta,item.par=data.frame(c(x$a),c(x$d)),BID=BID,Y=Y,I=SE,model=model,h=h)$sinI_vector
  I=matrix-tcrossprod(vector,vector)/nrow(theta)

  se <- matrix(sqrt(diag(solve(I))),nrow=max(BID$Block),ncol=5,byrow = TRUE)
  se.a <- se[,1:3]
  se.d <- se[,4:5]
  se.d <- cbind(se.d,sqrt(se.d[,1]^2+se.d[,2]^2))
  list(SE.a=se.a,SE.d=se.d,SE.sigm=NULL)
}


