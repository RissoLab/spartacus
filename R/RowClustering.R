RowClustering <- function(x, Ds, coordinates, iNN, n.neighbors, Alpha, Tau, Mu, Delta, Beta, Phi){
  K <- ifelse(is.vector(Mu), 1, nrow(Mu))
  goodK <- 1:K
  goodR <- unique(Ds)

  ll <- matrix(0, nrow(x), K)
  for (r in goodR) {
    for (k in goodK) {
      Linv <- computeLinv(covparms = c(1, Phi[r], Delta[k,r]),
                          covfun_name = "exponential_isotropic",
                          locs.ordered = coordinates[Ds == r, ],
                          iNN = iNN[[r]],
                          m = n.neighbors)

      nlogDet <- sum(log(diag(Linv)^2))
      Eta <- Linv%*%(t(x[, Ds == r])-Mu[k,r])
      Qv <- colSums(Eta*Eta)

      ll[,k] <- ll[,k] + (
        sum(Ds == r)*(log(Alpha-1)-log(Beta[k,r])-log(Tau))+
          nlogDet-(Alpha-1)*Qv/(Beta[k,r]*Tau)
      )
    }
  }
  allocation <- apply(ll,1,which.max)
  #stoch.allocation <- apply(ll,1,function(p){
  #    pp <- p-mean(p)
  #    sample(1:K, size = 1, prob = exp(pp))})
  return(list(allocation = allocation,
              #stoch.allocation = stoch.allocation,
              prob = ll))
}
