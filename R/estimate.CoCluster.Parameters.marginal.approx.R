estimate.CoCluster.Parameters.marginal.approx <- function(x,
                                                       locs.ordered,
                                                       iNN,
                                                       Alpha,
                                                       Tau,
                                                       mu0,
                                                       delta0,
                                                       beta0,
                                                       Phi
)
{
  n <- nrow(x)
  p <- ncol(x)
  cur.mu <- mu0
  cur.beta <- beta0
  cur.delta <- delta0

  # --update delta
  routine.delta <- optim(par = cur.delta, method = "Brent", fn = function(delta){
    Linv <- computeLinv(covparms = c(1, Phi, delta),
                        covfun_name = "exponential_isotropic",
                        locs.ordered = locs.ordered,
                        iNN = iNN)

    nlogDet <- 2*sum(log(diag(Linv)))
    Vinv <- t(Linv)%*%Linv
    mu.tmp <- as.numeric(colSums(Vinv)%*%colMeans(x)/sum(Vinv))
    Eta <- Linv%*%(t(x)-mu.tmp)
    Qt <- sum(Eta*Eta)
    return(n*p*log(Qt/(n*p))-n*nlogDet)
  }, lower = 1e-4, upper = 10^4, control = list(maxit = 20))
  if(routine.delta$convergence != 0){
    stop("Convergence error in delta!")
  }
  cur.delta <- routine.delta$par

  # -- update mu
  Linv <- computeLinv(covparms = c(1, Phi, cur.delta),
                      covfun_name = "exponential_isotropic",
                      locs.ordered = locs.ordered,
                      iNN = iNN)
  Vinv <- t(Linv)%*%Linv
  cur.mu <- as.numeric(colSums(Vinv)%*%colMeans(x)/sum(Vinv))

  # -- update beta
  Eta <- Linv%*%(t(x)-cur.mu)
  Qt <- sum(Eta*Eta)
  cur.beta = (Alpha-1)*Qt/(Tau*n*p)

  return(list(mu = cur.mu,
              delta = cur.delta,
              beta = cur.beta))
}
