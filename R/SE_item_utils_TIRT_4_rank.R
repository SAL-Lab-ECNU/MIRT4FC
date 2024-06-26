px_TIRT_4_rank <- function(theta,a,d,x=c(1,2,3,4)){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)
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
  if(identical(x,c(1,2,3,4))){       #顺序为1,2,3时的一阶导，下方同理
    L=pabcd/pall
  }else if(identical(x,c(1,2,4,3))){
    L=pabdc/pall
  }else if(identical(x,c(1,3,2,4))){
    L=pacbd/pall
  }else if(identical(x,c(1,3,4,2))){
    L=pacdb/pall
  }else if(identical(x,c(1,4,2,3))){
    L=padbc/pall
  }else if(identical(x,c(1,4,3,2))){
    L=padcb/pall
  }else if(identical(x,c(2,1,3,4))){
    L=pbacd/pall
  }else if(identical(x,c(2,1,4,3))){
    L=pbadc/pall
  }else if(identical(x,c(2,3,1,4))){
    L=pbcad/pall
  }else if(identical(x,c(2,3,4,1))){
    L=pbcda/pall
  }else if(identical(x,c(2,4,1,3))){
    L=pbdac/pall
  }else if(identical(x,c(2,4,3,1))){
    L=pbdca/pall
  }else if(identical(x,c(3,1,2,4))){
    L=pcabd/pall
  }else if(identical(x,c(3,1,4,2))){
    L=pcadb/pall
  }else if(identical(x,c(3,2,1,4))){
    L=pcbad/pall
  }else if(identical(x,c(3,2,4,1))){
    L=pcbda/pall
  }else if(identical(x,c(3,4,1,2))){
    L=pcdab/pall
  }else if(identical(x,c(3,4,2,1))){
    L=pcdba/pall
  }else if(identical(x,c(4,1,2,3))){
    L=pdabc/pall
  }else if(identical(x,c(4,1,3,2))){
    L=pdacb/pall
  }else if(identical(x,c(4,2,1,3))){
    L=pdbac/pall
  }else if(identical(x,c(4,2,3,1))){
    L=pdbca/pall
  }else if(identical(x,c(4,3,1,2))){
    L=pdcab/pall
  }else if(identical(x,c(4,3,2,1))){
    L=pdcba/pall
  }

  L
}




singleI_TIRT_4_rank <- function(theta,item.par,BID,Y,model,I,h=NULL){
  grad_FDM_item_TIRT_4_rank <- function(theta,a,d,x=c(1,2,3,4),h){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)

    g.ia <- (px_TIRT_4_rank(theta,a+c(h,0,0,0),d,x)-px_TIRT_4_rank(theta,a,d,x))/h
    g.ja <- (px_TIRT_4_rank(theta,a+c(0,h,0,0),d,x)-px_TIRT_4_rank(theta,a,d,x))/h
    g.ka <- (px_TIRT_4_rank(theta,a+c(0,0,h,0),d,x)-px_TIRT_4_rank(theta,a,d,x))/h
    g.la <- (px_TIRT_4_rank(theta,a+c(0,0,0,h),d,x)-px_TIRT_4_rank(theta,a,d,x))/h
    g.id <- (px_TIRT_4_rank(theta,a,d+c(h,0,0,-h),x)-px_TIRT_4_rank(theta,a,d,x))/h
    g.jd <- (px_TIRT_4_rank(theta,a,d+c(0,h,0,-h),x)-px_TIRT_4_rank(theta,a,d,x))/h
    g.kd <- (px_TIRT_4_rank(theta,a,d+c(0,0,h,-h),x)-px_TIRT_4_rank(theta,a,d,x))/h
    c(g.ia,g.ja,g.ka,g.la,g.id,g.jd,g.kd)
  }
  grad_CDM_item_TIRT_4_rank <- function(theta,a,d,x=c(1,2,3,4),h){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)

    g.ia <- (px_TIRT_4_rank(theta,a+c(h,0,0,0),d,x)-px_TIRT_4_rank(theta,a-c(h,0,0,0),d,x))/2*h
    g.ja <- (px_TIRT_4_rank(theta,a+c(0,h,0,0),d,x)-px_TIRT_4_rank(theta,a-c(0,h,0,0),d,x))/2*h
    g.ka <- (px_TIRT_4_rank(theta,a+c(0,0,h,0),d,x)-px_TIRT_4_rank(theta,a-c(0,0,h,0),d,x))/2*h
    g.la <- (px_TIRT_4_rank(theta,a+c(0,0,0,h),d,x)-px_TIRT_4_rank(theta,a-c(0,0,0,h),d,x))/2*h
    g.id <- (px_TIRT_4_rank(theta,a,d+c(h,0,0,-h),x)-px_TIRT_4_rank(theta,a,d-c(h,0,0,-h),x))/2*h
    g.jd <- (px_TIRT_4_rank(theta,a,d+c(0,h,0,-h),x)-px_TIRT_4_rank(theta,a,d-c(0,h,0,-h),x))/2*h
    g.kd <- (px_TIRT_4_rank(theta,a,d+c(0,0,h,-h),x)-px_TIRT_4_rank(theta,a,d-c(0,0,h,-h),x))/2*h
    c(g.ia,g.ja,g.ka,g.la,g.id,g.jd,g.kd)
  }
  grad_RES_item_TIRT_4_rank <- function(theta,a,d,x=c(1,2,3,4),h){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)

    g.ia <- (px_TIRT_4_rank(theta,a-c(2*h,0,0,0),d,x)+8*px_TIRT_4_rank(theta,a+c(h,0,0,0),d,x)
             -px_TIRT_4_rank(theta,a+c(2*h,0,0,0),d,x)-8*px_TIRT_4_rank(theta,a-c(h,0,0,0),d,x))/12*h
    g.ja <- (px_TIRT_4_rank(theta,a-c(0,2*h,0,0),d,x)+8*px_TIRT_4_rank(theta,a+c(0,h,0,0),d,x)
             -px_TIRT_4_rank(theta,a+c(0,2*h,0,0),d,x)-8*px_TIRT_4_rank(theta,a-c(0,h,0,0),d,x))/12*h
    g.ka <- (px_TIRT_4_rank(theta,a-c(0,0,2*h,0),d,x)+8*px_TIRT_4_rank(theta,a+c(0,0,h,0),d,x)
             -px_TIRT_4_rank(theta,a+c(0,0,2*h,0),d,x)-8*px_TIRT_4_rank(theta,a-c(0,0,h,0),d,x))/12*h
    g.la <- (px_TIRT_4_rank(theta,a-c(0,0,0,2*h),d,x)+8*px_TIRT_4_rank(theta,a+c(0,0,0,h),d,x)
             -px_TIRT_4_rank(theta,a+c(0,0,0,2*h),d,x)-8*px_TIRT_4_rank(theta,a-c(0,0,0,h),d,x))/12*h
    g.id <- (px_TIRT_4_rank(theta,a,d-c(2*h,0,0,-2*h),x)+8*px_TIRT_4_rank(theta,a,d+c(h,0,0,-h),x)
             -px_TIRT_4_rank(theta,a,d+c(2*h,0,0,-2*h),x)-8*px_TIRT_4_rank(theta,a,d-c(h,0,0,-h),x))/12*h
    g.jd <- (px_TIRT_4_rank(theta,a,d-c(0,2*h,0,-2*h),x)+8*px_TIRT_4_rank(theta,a,d+c(0,h,0,-h),x)
             -px_TIRT_4_rank(theta,a,d+c(0,2*h,0,-2*h),x)-8*px_TIRT_4_rank(theta,a,d-c(0,h,0,-h),x))/12*h
    g.kd <- (px_TIRT_4_rank(theta,a,d-c(0,0,2*h,-2*h),x)+8*px_TIRT_4_rank(theta,a,d+c(0,0,h,-h),x)
             -px_TIRT_4_rank(theta,a,d+c(0,0,2*h,-2*h),x)-8*px_TIRT_4_rank(theta,a,d-c(0,0,h,-h),x))/12*h
    c(g.ia,g.ja,g.ka,g.la,g.id,g.jd,g.kd)
  }





  FDMEss_item_TIRT_4_rank <- function(theta,a,d,Y,h=1e-5){
    all.x <- list(c(1,2,3,4),c(1,2,4,3),c(1,3,2,4),c(1,3,4,2),c(1,4,2,3),c(1,4,3,2),
                  c(2,1,3,4),c(2,1,4,3),c(2,3,1,4),c(2,3,4,1),c(2,4,1,3),c(2,4,3,1),
                  c(3,1,2,4),c(3,1,4,2),c(3,2,1,4),c(3,2,4,1),c(3,4,1,2),c(3,4,2,1),
                  c(4,1,2,3),c(4,1,3,2),c(4,2,1,3),c(4,2,3,1),c(4,3,1,2),c(4,3,2,1))
    g <- grad_FDM_item_TIRT_4_rank(theta,a,d,x=all.x[[Y]],h=h)
    ret <- tcrossprod(g,g)/px_TIRT_4_rank(theta,a,d,xx=all.x[[Y]])^2

    ret
  }
  RESEss_item_TIRT_4_rank <- function(theta,a,d,Y,h=1){
    all.x <- list(c(1,2,3,4),c(1,2,4,3),c(1,3,2,4),c(1,3,4,2),c(1,4,2,3),c(1,4,3,2),
                  c(2,1,3,4),c(2,1,4,3),c(2,3,1,4),c(2,3,4,1),c(2,4,1,3),c(2,4,3,1),
                  c(3,1,2,4),c(3,1,4,2),c(3,2,1,4),c(3,2,4,1),c(3,4,1,2),c(3,4,2,1),
                  c(4,1,2,3),c(4,1,3,2),c(4,2,1,3),c(4,2,3,1),c(4,3,1,2),c(4,3,2,1))
    g <- grad_RES_item_TIRT_4_rank(theta,a,d,x=all.x[[Y]],h=h)
    ret <- tcrossprod(g,g)/px_TIRT_4_rank(theta,a,d,xx=all.x[[Y]])^2

    ret
  }
  CDMEss_item_TIRT_4_rank <- function(theta,a,d,Y,h=1){
    all.x <- list(c(1,2,3,4),c(1,2,4,3),c(1,3,2,4),c(1,3,4,2),c(1,4,2,3),c(1,4,3,2),
                  c(2,1,3,4),c(2,1,4,3),c(2,3,1,4),c(2,3,4,1),c(2,4,1,3),c(2,4,3,1),
                  c(3,1,2,4),c(3,1,4,2),c(3,2,1,4),c(3,2,4,1),c(3,4,1,2),c(3,4,2,1),
                  c(4,1,2,3),c(4,1,3,2),c(4,2,1,3),c(4,2,3,1),c(4,3,1,2),c(4,3,2,1))
    g <- grad_CDM_item_TIRT_4_rank(theta,a,d,x=all.x[[Y]],h=h)
    ret <- tcrossprod(g,g)/px_TIRT_4_rank(theta,a,d,xx=all.x[[Y]])^2

    ret
  }



  FDMSS_item_TIRT_4_rank <- function(theta,a,d,Y,h=1e-5){
    all.x <- list(c(1,2,3,4),c(1,2,4,3),c(1,3,2,4),c(1,3,4,2),c(1,4,2,3),c(1,4,3,2),
                  c(2,1,3,4),c(2,1,4,3),c(2,3,1,4),c(2,3,4,1),c(2,4,1,3),c(2,4,3,1),
                  c(3,1,2,4),c(3,1,4,2),c(3,2,1,4),c(3,2,4,1),c(3,4,1,2),c(3,4,2,1),
                  c(4,1,2,3),c(4,1,3,2),c(4,2,1,3),c(4,2,3,1),c(4,3,1,2),c(4,3,2,1))
    g <- grad_FDM_item_TIRT_4_rank(theta,a,d,x=all.x[[Y]],h=h)
    g
  }
  RESSS_item_TIRT_4_rank <- function(theta,a,d,Y,h=1){
    all.x <- list(c(1,2,3,4),c(1,2,4,3),c(1,3,2,4),c(1,3,4,2),c(1,4,2,3),c(1,4,3,2),
                  c(2,1,3,4),c(2,1,4,3),c(2,3,1,4),c(2,3,4,1),c(2,4,1,3),c(2,4,3,1),
                  c(3,1,2,4),c(3,1,4,2),c(3,2,1,4),c(3,2,4,1),c(3,4,1,2),c(3,4,2,1),
                  c(4,1,2,3),c(4,1,3,2),c(4,2,1,3),c(4,2,3,1),c(4,3,1,2),c(4,3,2,1))
    g <- grad_RES_item_TIRT_4_rank(theta,a,d,x=all.x[[Y]],h=h)
    g
  }
  CDMSS_item_TIRT_4_rank <- function(theta,a,d,Y,h=1){
    all.x <- list(c(1,2,3,4),c(1,2,4,3),c(1,3,2,4),c(1,3,4,2),c(1,4,2,3),c(1,4,3,2),
                  c(2,1,3,4),c(2,1,4,3),c(2,3,1,4),c(2,3,4,1),c(2,4,1,3),c(2,4,3,1),
                  c(3,1,2,4),c(3,1,4,2),c(3,2,1,4),c(3,2,4,1),c(3,4,1,2),c(3,4,2,1),
                  c(4,1,2,3),c(4,1,3,2),c(4,2,1,3),c(4,2,3,1),c(4,3,1,2),c(4,3,2,1))
    g <- grad_CDM_item_TIRT_4_rank(theta,a,d,x=all.x[[Y]],h=h)
    g
  }


  itemI_TIRT_4_rank <- function(theta,item.par,BID,Y,I,h=NULL){
    D <- max(BID$Block)*(max(BID$Item)*2-1) #number of dims
    itemIEss <- diag(rep(0,D))
    itemISS <- rep(0,D)
    for(b in unique(BID$Block)){
      items <- which(BID$Block==b)
      item <- c((b-1)*7+1,(b-1)*7+2,(b-1)*7+3,(b-1)*7+4,(b-1)*7+5,(b-1)*7+6,(b-1)*7+7)
      dim.b <- BID$Dim[items]
      a.b <- item.par[items,1]
      d.b <- item.par[items,2]
      theta.b <- theta[dim.b]
      Y.b <- Y[b]
      if(I=="XPD"){
        itemIEss[item,item] <- FDMEss_item_TIRT_4_rank(theta.b,a.b,d.b,Y.b)
        itemISS[item] <- FDMSS_item_TIRT_4_rank(theta.b,a.b,d.b,Y.b)
      }else if(I=="FDM"){
        if(is.null(h)){
          itemIEss[item,item] <- FDMEss_item_TIRT_4_rank(theta.b,a.b,d.b,Y.b)
          itemISS[item] <- FDMSS_item_TIRT_4_rank(theta.b,a.b,d.b,Y.b)
        }else{
          itemIEss[item,item] <- FDMEss_item_TIRT_4_rank(theta.b,a.b,d.b,Y.b,h=h)
          itemISS[item] <- FDMSS_item_TIRT_4_rank(theta.b,a.b,d.b,Y.b,h=h)
        }
      }else if(I=="CDM"){
        if(is.null(h)){
          itemIEss[item,item] <- CDMEss_item_TIRT_4_rank(theta.b,a.b,d.b,Y.b)
          itemISS[item] <- CDMSS_item_TIRT_4_rank(theta.b,a.b,d.b,Y.b)
        }else{
          itemIEss[item,item] <- CDMEss_item_TIRT_4_rank(theta.b,a.b,d.b,Y.b,h=h)
          itemISS[item] <- CDMSS_item_TIRT_4_rank(theta.b,a.b,d.b,Y.b,h=h)
        }
      }else if(I=="RES"){
        if(is.null(h)){
          itemIEss[item,item] <- RESEss_item_TIRT_4_rank(theta.b,a.b,d.b,Y.b)
          itemISS[item] <- RESSS_item_TIRT_4_rank(theta.b,a.b,d.b,Y.b)
        }else{
          itemIEss[item,item] <- RESEss_item_TIRT_4_rank(theta.b,a.b,d.b,Y.b,h=h)
          itemISS[item] <- RESSS_item_TIRT_4_rank(theta.b,a.b,d.b,Y.b,h=h)
        }
      }
    }
    list(itemIEss=itemIEss,itemISS=itemISS)
  }


  # information matrix of block b

  ## 入参：以block为单位计算
  ## theta=作答者能力值，a=题目的a参数，d=（ a * b）





  sinI_matrix=foreach(i = seq_len(nrow(theta)),.inorder=FALSE,.combine='+',.export = c("foreach","%dopar%","px_TIRT_4_rank"))%dopar%{
                                                                                         itemI_TIRT_4_rank(theta[i,],item.par,BID,Y[i,],I,h=h)$itemIEss
                                                                                       }
  sinI_vector=foreach(i = seq_len(nrow(theta)),.inorder=FALSE,.combine='+',.export = c("foreach","%dopar%","px_TIRT_4_rank"))%dopar%{
                                                                                         itemI_TIRT_4_rank(theta[i,],item.par,BID,Y[i,],I,h=h)$itemISS
                                                                                       }
  list(sinI_matrix=sinI_matrix,sinI_vector=sinI_vector)
}

itemSE_TIRT_4_rank <- function(x,BID,Y,model,SE,IN,theta,h=NULL){


  matrix=singleI_TIRT_4_rank(theta=theta,item.par=data.frame(c(x$a),c(x$d)),BID=BID,Y=Y,I=SE,model=model,h=h)$sinI_matrix
  vector=singleI_TIRT_4_rank(theta=theta,item.par=data.frame(c(x$a),c(x$d)),BID=BID,Y=Y,I=SE,model=model,h=h)$sinI_vector
  I=matrix-tcrossprod(vector,vector)/nrow(theta)

  se <- matrix(sqrt(diag(solve(I))),nrow=max(BID$Block),ncol=7,byrow = TRUE)
  se.a <- se[,1:4]
  se.d <- se[,5:7]
  se.d <- cbind(se.d,sqrt(se.d[,1]^2+se.d[,2]^2+se.d[,3]^2))
  list(SE.a=se.a,SE.d=se.d,SE.sigm=NULL)
}


