data.sim_ZG <- function(item.par=item.par,theta=theta,BID=BID,blocksize=3,res='rank'){
  if(blocksize==2){
    Y <- data.sim_ZG_2(item.par=item.par,theta=theta,BID=BID)
  }else{
    cat("The model is not included")
  }
  Y
}
