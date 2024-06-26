theta.est_TIRT <- function(Y,a=NULL,d=NULL,BID,sigma=NULL,prior=TRUE,blocksize=3,res='rank',model='TIRT'){
  D <- max(BID$Dim)
  N <- nrow(Y)
  theta.est <- matrix(0,N,D)
  theta.se <- matrix(NA,N,D)
  if(blocksize==2){

    for (i in 1:N){
      theta.est[i,] <- lbfgs(rep(0,D),fn=logL_TIRT_2,lower = rep(-5,D),upper = rep(5,D),a=a,d=d,yi=Y[i,],BID=BID,sigma=sigma)$par
      #theta.est[i,] <- optim(rep(0,D),fn = logLi_TIRT_2,method = "L-BFGS-B",lower = rep(-5,D), upper = rep(5,D),a=a, d=d,yi=Y[i,], BID=BID,prior=prior,sigma=sigma,control=list(fnscale=-1))$par
      theta.se[i,] <- test.info_theta_TIRT_2(theta.est[i,], item.par = data.frame(a=c(a),d=c(d)), BID = BID)$se
    }

  }else if(blocksize==3){
    if(res=='rank'|res=='mole'){

      for (i in 1:N){
        theta.est[i,] <- lbfgs(rep(0,D),fn=logL_TIRT_3_rank,lower = rep(-5,D),upper = rep(5,D),a=a,d=d,yi=Y[i,],BID=BID,sigma=sigma)$par
        #theta.est[i,] <- optim(rep(0,D),fn = logLi_TIRT_3_rank,method = "L-BFGS-B",lower = rep(-5,D), upper = rep(5,D),a=a, d=d,yi=Y[i,], BID=BID,prior=prior,sigma=sigma,control=list(fnscale=-1))$par
        theta.se[i,] <- test.info_theta_TIRT_3_rank(theta.est[i,], item.par = data.frame(a=c(a),d=c(d)), BID = BID)$se
      }


    }else if(res=='pick'){

      for (i in 1:N){
        theta.est[i,] <- lbfgs(rep(0,D),fn=logL_TIRT_3_pick,lower = rep(-5,D),upper = rep(5,D),a=a,d=d,yi=Y[i,],BID=BID,sigma=sigma)$par
        #theta.est[i,] <- optim(rep(0,D),fn = logLi_TIRT_3_pick,method = "L-BFGS-B",lower = rep(-5,D), upper = rep(5,D),a=a, d=d,yi=Y[i,], BID=BID,prior=prior,sigma=sigma,control=list(fnscale=-1))$par
        theta.se[i,] <- test.info_theta_TIRT_3_pick(theta.est[i,], item.par = data.frame(a=c(a),d=c(d)), BID = BID)$se
      }
    }
  }else if(blocksize==4){
    if(res=='rank'){

      for (i in 1:N){
        theta.est[i,] <- lbfgs(rep(0,D),fn=logL_TIRT_4_rank,lower = rep(-5,D),upper = rep(5,D),a=a,d=d,yi=Y[i,],BID=BID,sigma=sigma)$par
        #theta.est[i,] <- optim(rep(0,D),fn = logLi_TIRT_4_rank,method = "L-BFGS-B",lower = rep(-5,D), upper = rep(5,D),a=a, d=d,yi=Y[i,], BID=BID,prior=prior,sigma=sigma,control=list(fnscale=-1))$par
        theta.se[i,] <- test.info_theta_TIRT_4_rank(theta.est[i,], item.par = data.frame(a=c(a),d=c(d)), BID = BID)$se
      }
    }else if(res=='pick'){

      for (i in 1:N){
        theta.est[i,] <- lbfgs(rep(0,D),fn=logL_TIRT_4_pick,lower = rep(-5,D),upper = rep(5,D),a=a,d=d,yi=Y[i,],BID=BID,sigma=sigma)$par
        #theta.est[i,] <- optim(rep(0,D),fn = logLi_TIRT_4_pick,method = "L-BFGS-B",lower = rep(-5,D), upper = rep(5,D),a=a, d=d,yi=Y[i,], BID=BID,prior=prior,sigma=sigma,control=list(fnscale=-1))$par
        theta.se[i,] <- test.info_theta_TIRT_4_pick(theta.est[i,], item.par = data.frame(a=c(a),d=c(d)), BID = BID)$se
      }
    }else if(res=='mole'){

      for (i in 1:N){
        theta.est[i,] <- lbfgs(rep(0,D),fn=logL_TIRT_4_mole,lower = rep(-5,D),upper = rep(5,D),a=a,d=d,yi=Y[i,],BID=BID,sigma=sigma)$par
        #theta.est[i,] <- optim(rep(0,D),fn = logLi_TIRT_4_mole,method = "L-BFGS-B",lower = rep(-5,D), upper = rep(5,D),a=a, d=d,yi=Y[i,], BID=BID,prior=prior,sigma=sigma,control=list(fnscale=-1))$par
        theta.se[i,] <- test.info_theta_TIRT_4_mole(theta.est[i,], item.par = data.frame(a=c(a),d=c(d)), BID = BID)$se
      }
    }
  }



  return(list(theta.est=theta.est,theta.se=theta.se))
}




