
px_item_2PL_2 <- function(theta,a,d,x=1){#一个block内的一阶导，#a=c(ai,aj,ak) d=c(di,dj,dk) theta=c(thetai,thetaj,thetak)
  z <- exp(t(a*t(theta)-d))          #1*3

  # p.j.ijk <- p.j.ikj <- p.j.jki <- p.j.jik <- p.j.kij <- p.j.kji
  p.j <- z[,2]/rowSums(z)
  # p.i.ijk <- p.i.ikj <- p.i.jki <- p.i.jik <- p.i.kij <- p.i.kji
  p.i <- z[,1]/rowSums(z)

  if(identical(x,1)){       #顺序为1,2,3时的一阶导，下方同理
    L=p.i
  }else if(identical(x,2)){
    L=p.j
  }
  L
}



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
    z <- exp(t(a*t(theta)-d))
    thetai <- theta[,1]
    thetaj <- theta[,2]

    # p.j.ijk <- p.j.ikj <- p.j.jki <- p.j.jik <- p.j.kij <- p.j.kji
    p.j <- z[,2]/rowSums(z)
    # p.i.ijk <- p.i.ikj <- p.i.jki <- p.i.jik <- p.i.kij <- p.i.kji
    p.i <- z[,1]/rowSums(z)

    if(identical(x,1)){
      g.ia <- thetai*p.i*p.j
      g.ja <- -thetaj*p.i*p.j
      g.id <- -2*p.i*p.j
    }else if(identical(x,2)){
      g.ja <- thetaj*p.i*p.j
      g.ia <- -thetai*p.j*p.i
      g.id <- 2*p.i*p.j
    }

    cbind(g.ia,g.ja,g.id)
  }



  FDM_item_2PL_2 <- function(theta,a,d,Y,h=1e-5){
    all.x <- list(1,2)
    g <- grad_FDM_item_2PL_2(theta,a,d,x=all.x[[Y]],h=h)/px_2PL_2(theta,a,d,x=all.x[[Y]])
    g
  }
  RES_item_2PL_2 <- function(theta,a,d,Y,h=1){
    all.x <- list(1,2)
    g <- grad_RES_item_2PL_2(theta,a,d,x=all.x[[Y]],h=h)/px_2PL_2(theta,a,d,x=all.x[[Y]])
    g
  }
  CDM_item_2PL_2 <- function(theta,a,d,Y,h=1){
    all.x <- list(1,2)
    g <- grad_CDM_item_2PL_2(theta,a,d,x=all.x[[Y]],h=h)/px_2PL_2(theta,a,d,x=all.x[[Y]])
    g
  }
  XPD_item_2PL_2 <- function(theta,a,d,Y){
    all.x <- list(1,2)
    g <- grad_item_2PL_2(theta,a,d,x=all.x[[Y]])/px_item_2PL_2(theta,a,d,x=all.x[[Y]])
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
    itemI <- matrix(0,nrow=nrow(theta),ncol=D)
    for(b in unique(BID$Block)){
      items <- which(BID$Block==b)
      item <- c((b-1)*3+1,(b-1)*3+2,(b-1)*3+3)
      dim.b <- BID$Dim[items]
      a.b <- item.par[items,1]
      d.b <- item.par[items,2]
      theta.b <- theta[,dim.b]
      Y.b <- Y[b]
      if(I=="XPD"){
        itemI[,item] <- XPD_item_2PL_2(theta.b,a.b,d.b,Y.b)
      }else if(I=="FDM"){
        if(is.null(h)){
          itemI[item] <- FDM_item_2PL_2(theta.b,a.b,d.b,Y.b)
        }else{
          itemI[item] <- FDM_item_2PL_2(theta.b,a.b,d.b,Y.b,h=h)
        }
      }else if(I=="CDM"){
        if(is.null(h)){
          itemI[item] <- CDM_item_2PL_2(theta.b,a.b,d.b,Y.b)
        }else{
          itemI[item] <- CDM_item_2PL_2(theta.b,a.b,d.b,Y.b,h=h)
        }
      }else if(I=="RES"){
        if(is.null(h)){
          itemI[item] <- RES_item_2PL_2(theta.b,a.b,d.b,Y.b)
        }else{
          itemI[item] <- RES_item_2PL_2(theta.b,a.b,d.b,Y.b,h=h)
        }
      }
    }
    itemI
  }


  # information matrix of block b

  ## 入参：以block为单位计算
  ## theta=作答者能力值，a=题目的a参数，d=（ a * b）
  itemSE_Fisher_2PL_2 <- function(x,BID,Y,IN,q=10,h=NULL){
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

    for (h in 1:nrow(resorder)) {
      newresorder1=matrix(rep(as.matrix(resorder)[h,], each = q^D), nrow = q^D, ncol = J, byrow = F)
      for(i in 1:q^D){
        for(j in 1:J){
          p1_5[i,j]=px_2PL_2(X_matrix[i,],a[,j],d[,j],x=as.numeric(newresorder1[i,j]))
        }
      }
      paiH1[h,h]=sum(apply(p1_5, 1, prod)*AX_p)
      post=apply(p1_5, 1, prod)*AX_p/paiH1[h,h]

      sinI_vector[,h]=foreach(i = seq_len(nrow(X_matrix)),.inorder=FALSE,.combine='+',.export = c("itemI_2PL_2",                                                                                                                                                                                    "grad_item_2PL_2",
                                                                                                  "px_2PL_2",
                                                                                                  "XPD_item_2PL_2"))%dopar%{
                                                                                                    itemI_2PL_2(X_matrix[i,],data.frame(a=c(x$a),d=c(x$d)),BID,resorder[h,],I='XPD',h=NULL)*post[i]
                                                                                                  }
    }
    I=sinI_vector%*%paiH1%*%t(sinI_vector)*N
    se.a=NULL
    se.d=NULL
    se.sigm=NULL
    tryCatch({
      se <- matrix(sqrt(diag(solve(I))),nrow=max(BID$Block),ncol=IN,byrow = TRUE)
      se.a <- se[,1:2]
      se.d <- se[,3]
      se.d <- cbind(se.d,se.d)
    }, error = function(e) {
      cat("The information matrix is irreversible!")
    })
    list(SE.a=se.a,SE.d=se.d,SE.sigm=se.sigm,I=I)
  }



  itemSE_Simple_2PL_2 <- function(x,BID,Y,SE='XPD',IN,hess=NULL,h1=NULL){
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
    for (h in 1:nrow(Y)) {
      if(h %% 200==0){
        cat(h,'/',nrow(Y),'\n')
      }
      step_value =nrow(x$thetalist[[1]])
      step_length=length(x$thetalist)
      # 生成新向量
      new_vector <- sapply(which(apply(Y, 1, function(x) all(x == as.matrix(Y[h,])))), function(x) seq(from = x, by = step_value, length.out = step_length))
      thetaexchoose=thetaex[c(new_vector),]


      if(SE=='Louis'|SE=='Sandwich'){
        g=itemI_2PL_2(thetaexchoose,data.frame(a=c(x$a),d=c(x$d)),BID,as.matrix(Y[h,]),I='XPD',h=NULL)
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
        sinI_vector[,h]=colMeans(itemI_2PL_2(thetaexchoose,data.frame(a=c(x$a),d=c(x$d)),BID,as.matrix(Y[h,]),I=SE,h=h1))

      }

    }

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
      se.a <- se[,1:2]
      se.d <- se[,3]
      se.d <- cbind(se.d,se.d)
    }, error = function(e) {
      cat("The information matrix is irreversible!")
    })
    s2=Sys.time()
    list(SE.a=se.a,SE.d=se.d,SE.sigm=se.sigm,I=I,time=s2-s1)
  }


