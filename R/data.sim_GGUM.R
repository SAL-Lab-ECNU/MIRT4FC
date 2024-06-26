data.sim_GGUM <- function(item.par=item.par,theta=theta,BID=BID,blocksize=3,res='rank'){
  if(blocksize==2){
    Y <- data.sim_GGUM_2(item.par=item.par,theta=theta,BID=BID)
  }else if(blocksize==3){
    if(res=='pick'){
      Y <- data.sim_GGUM_3_pick(item.par=item.par,theta=theta,BID=BID)
    }else if(res=='rank'|res=='mole'){
      Y <- data.sim_GGUM_3_rank(item.par=item.par,theta=theta,BID=BID)
    }
  }else if(blocksize==4){
    if(res=='pick'){
      Y <- data.sim_GGUM_4_pick(item.par=item.par,theta=theta,BID=BID)
    }else if(res=='rank'){
      Y <- data.sim_GGUM_4_rank(item.par=item.par,theta=theta,BID=BID)
    }else if(res=='mole'){
      Y <- data.sim_GGUM_4_mole(item.par=item.par,theta=theta,BID=BID)
    }
  }
  Y
}
