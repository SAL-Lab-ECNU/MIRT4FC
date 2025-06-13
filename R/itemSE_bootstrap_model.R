itemSE_bootstrap_2PL<- function(plist,blocksize,IN,J,D){

  SE.a=NULL
  SE.d=NULL
  SE.sigm=NULL
  SE_big <- NULL
  for (i in 1:length(plist)) {
    if (is.null(SE_big)) {
      SE_big <- plist[[i]]
    } else {
      SE_big <- cbind(SE_big, plist[[i]])
    }
  }
  SE.all=apply(SE_big, 1, sd)
  for (j in 1:J) {
    start_index <- (j - 1) * IN + 1
    end_index <- min(start_index + blocksize - 1, length(SE.all))
    extracted_part <- SE.all[start_index:end_index]
    SE.a <- c(SE.a, extracted_part)
  }
  SE.d_miss=NULL
  for (j in 1:J) {
    start_index <- (j - 1) * IN + 1+blocksize
    end_index <- min(start_index + blocksize - 2, length(SE.all))
    extracted_part <- SE.all[start_index:end_index]
    SE.d_miss <- c(SE.d_miss, extracted_part)
  }

  for (i in seq(1, length(SE.d_miss), by = blocksize-1)) {
    if(blocksize==2){
      third_para <- sqrt(SE.d_miss[i]^2)
      SE.d <- c(SE.d, SE.d_miss[i], third_para)
    }else if(blocksize==3){
      third_para <- sqrt(SE.d_miss[i]^2+SE.d_miss[i+1]^2)
      SE.d <- c(SE.d, SE.d_miss[i], SE.d_miss[i+1], third_para)
    }else if(blocksize==4){
      third_para <- sqrt(SE.d_miss[i]^2+SE.d_miss[i+1]^2+SE.d_miss[i+2]^2)
      SE.d <- c(SE.d, SE.d_miss[i], SE.d_miss[i+1], SE.d_miss[i+2], third_para)
    }
  }
  #改变sigm的SE格式
  SE.sigm_miss=SE.all[(J*IN+1):length(SE.all)]
  SE.sigm <- matrix(0, nrow = D, ncol = D)
  index=1
  for (i in 2:D) {
    for (j in 1:(i-1)) {
      SE.sigm[i, j] <- SE.sigm_miss[index]
      index <- index+1
    }
  }
  SE.sigm[upper.tri(SE.sigm)] <- t(SE.sigm)[upper.tri(SE.sigm)]#填充上三角
  return(list(SE.a=matrix(SE.a,ncol=J,nrow=blocksize),SE.d=matrix(SE.d,ncol=J,nrow=blocksize),SE.sigm=SE.sigm))
}

