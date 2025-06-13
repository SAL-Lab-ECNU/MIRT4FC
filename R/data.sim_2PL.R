data.sim_2PL <- function(item.par=item.par,theta=theta,BID=BID,blocksize=3,res='rank'){
  if(blocksize==2){
    Y <- data.sim_2PL_2(item.par=item.par,theta=theta,BID=BID)
  }else if(blocksize==3){
    if(res=='pick'){
      Y <- data.sim_2PL_3_pick(item.par=item.par,theta=theta,BID=BID)
    }else if(res=='rank'|res=='mole'){
      Y <- data.sim_2PL_3_rank(item.par=item.par,theta=theta,BID=BID)
    }
  }else if(blocksize==4){
    if(res=='pick'){
      Y <- data.sim_2PL_4_pick(item.par=item.par,theta=theta,BID=BID)
    }else if(res=='rank'){
      Y <- data.sim_2PL_4_rank(item.par=item.par,theta=theta,BID=BID)
    }else if(res=='mole'){
      Y <- data.sim_2PL_4_mole(item.par=item.par,theta=theta,BID=BID)
    }
  }
  Y
}
