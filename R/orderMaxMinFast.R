#From NPVecchia
orderMaxMinFast <- function( locs, numpropose){

  n <- nrow(locs)
  d <- ncol(locs)
  remaininginds <- 1:n
  orderinds <- rep(0L,n)
  # pick a center point
  mp <- matrix(colMeans(locs),1,d)
  distmp <- fields::rdist(locs,mp)
  ordermp <- order(distmp)
  orderinds[1] = ordermp[1]
  remaininginds <- remaininginds[remaininginds!=orderinds[1]]
  for( j in 2:(n-1) ){
    randinds <- sample(remaininginds,min(numpropose,length(remaininginds)))
    distarray <-  fields::rdist(locs[orderinds[1:j-1],,drop=FALSE],locs[randinds,,drop=FALSE])
    bestind <- which(colMins(distarray, useNames = F) ==  max( colMins( distarray, useNames = F) ))
    orderinds[j] <- randinds[bestind[1]]
    remaininginds <- remaininginds[remaininginds!=orderinds[j]]
  }
  orderinds[n] <- remaininginds
  orderinds
}
