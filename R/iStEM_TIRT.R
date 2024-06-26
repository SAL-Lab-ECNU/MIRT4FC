iStEM_TIRT <- function(Y,BID,positive=rep(TRUE,nrow(BID)),blocksize=3,res='rank',M=10,B=20,SE='Bootstrap',
                      sigma=NULL,theta=NULL,fix.sigma = FALSE,burnin.maxitr=40,model='TIRT',
                      maxitr=500,eps1=1.5,eps2=0.4,frac1=.2,frac2=.5,cores=NULL,h=NULL){
  if(blocksize==2){
    x <- iStEM_TIRT_2(Y=Y,BID=BID,positive=positive,M=M,B=B,SE=SE,blocksize=blocksize,model=model,
                     sigma=sigma,theta=theta,fix.sigma = fix.sigma,burnin.maxitr=burnin.maxitr,
                     maxitr=maxitr,eps1=eps1,eps2=eps2,frac1=frac1,frac2=frac2,cores=cores,h=h)
  }else if(blocksize==3){
    if(res=='pick'){
      x <- iStEM_TIRT_3_pick(Y=Y,BID=BID,positive=positive,M=M,B=B,SE=SE,blocksize=blocksize,model=model,
                            sigma=sigma,theta=theta,fix.sigma = fix.sigma,burnin.maxitr=burnin.maxitr,
                            maxitr=maxitr,eps1=eps1,eps2=eps2,frac1=frac1,frac2=frac2,cores=cores,h=h)
    }else if(res=='rank'|res=='mole'){
      x <- iStEM_TIRT_3_rank(Y=Y,BID=BID,positive=positive,M=M,B=B,SE=SE,blocksize=blocksize,model=model,
                            sigma=sigma,theta=theta,fix.sigma = fix.sigma,burnin.maxitr=burnin.maxitr,
                            maxitr=maxitr,eps1=eps1,eps2=eps2,frac1=frac1,frac2=frac2,cores=cores,h=h)
    }
  }else if(blocksize==4){
    if(res=='pick'){
      x <- iStEM_TIRT_4_pick(Y=Y,BID=BID,positive=positive,M=M,B=B,SE=SE,blocksize=blocksize,model=model,
                            sigma=sigma,theta=theta,fix.sigma = fix.sigma,burnin.maxitr=burnin.maxitr,
                            maxitr=maxitr,eps1=eps1,eps2=eps2,frac1=frac1,frac2=frac2,cores=cores,h=h)
    }else if(res=='rank'){
      x <- iStEM_TIRT_4_rank(Y=Y,BID=BID,positive=positive,M=M,B=B,SE=SE,blocksize=blocksize,model=model,
                            sigma=sigma,theta=theta,fix.sigma = fix.sigma,burnin.maxitr=burnin.maxitr,
                            maxitr=maxitr,eps1=eps1,eps2=eps2,frac1=frac1,frac2=frac2,cores=cores,h=h)
    }else if(res=='mole'){
      x <- iStEM_TIRT_4_mole(Y=Y,BID=BID,positive=positive,M=M,B=B,SE=SE,blocksize=blocksize,model=model,
                            sigma=sigma,theta=theta,fix.sigma = fix.sigma,burnin.maxitr=burnin.maxitr,
                            maxitr=maxitr,eps1=eps1,eps2=eps2,frac1=frac1,frac2=frac2,cores=cores,h=h)
    }
  }
  return(x)
}
