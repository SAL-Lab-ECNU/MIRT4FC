batch.var <- function(plist,n){
  phi.bar <- sapply(plist,rowMeans) # of par x n
  phi.hat <- rowMeans(Reduce(cbind,plist))
  rowMeans((phi.bar-phi.hat)^2)/(n-1)
}
