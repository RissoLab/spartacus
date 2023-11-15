estimatePhi <- function(x, Cs, locs.ordered, iNN, Alpha, Tau, Mu, Delta, Beta, phi.old = 1, n.neighbors){
  goodK <- sort(unique(Cs))

  val <- rep(0,length(goodK))
  routine.phi <- optim(par = phi.old, method = "Brent", fn = function(phi){
    for(k in goodK){
      Linv <- computeLinv(covparms = c(1, phi, Delta[k]),
                          covfun_name = "exponential_isotropic",
                          locs.ordered = locs.ordered,
                          iNN = iNN,
                          m = n.neighbors)
      nlogDet <- 2*sum(log(diag(Linv)))
      Eta <- Linv%*%(t(x[Cs == k,])-Mu[k])
      Qt <- sum(Eta*Eta)
      val[k] <- ncol(x[Cs == k,])*nrow(x[Cs == k,])*(log(Alpha-1)-log(Beta[k])-log(Tau)) +
        nrow(x[Cs == k,])*nlogDet-(Alpha-1)*Qt/(Beta[k]*Tau)
    }
    return(-sum(val))
  }, lower = 1e-4, upper = 1000)
  if(routine.phi$convergence != 0){
    stop("Convergence error in phi!")
  }
  return(routine.phi$par)
}

