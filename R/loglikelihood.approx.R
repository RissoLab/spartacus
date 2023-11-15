logL.Cocluster.approx <- function(x, locs.ordered, iNN, n.neighbors, Alpha, Tau, Mu, Delta, Beta, Phi){
  if(is.vector(x)) x <- matrix(x, nrow = 1)

  Linv <- computeLinv(covparms = c(1, Phi, Delta),
                      covfun_name = "exponential_isotropic",
                      locs.ordered = locs.ordered,
                      iNN = iNN,
                      m = n.neighbors)
  nlogDet <- 2*sum(log(diag(Linv)))
  Eta <- Linv%*%(t(x)-Mu)
  Qt <- sum(Eta*Eta)

  #logL
  ncol(x)*nrow(x)*(log(Alpha-1)-log(Beta)-log(Tau)) +
    nrow(x)*nlogDet-(Alpha-1)*Qt/(Beta*Tau)
}
