updateAlphaBeta.CoCluster <- function(x,
                                      locs.ordered,
                                      Linv = NULL,
                                      iNN = NULL,
                                      n.neighbors = NULL,
                                      Alpha,
                                      Tau,
                                      Mu,
                                      Delta,
                                      Beta,
                                      Phi,
                                      max.iter = 1000,
                                      conv.criterion = list(tollerance = 1e-3, iterations = 10)
){
  n <- nrow(x)
  p <- ncol(x)
  conv.counter <- 0

  Quadratic <- rep(NA, n)
  if(is.null(Linv)){
    Linv <- computeLinv(covparms = c(1, Phi, Delta),
                        covfun_name = "exponential_isotropic",
                        locs.ordered = locs.ordered,
                        iNN = iNN,
                        m = n.neighbors)
  }
  nlogDet <- 2*sum(log(diag(Linv)))
  Eta <- Linv%*%(t(x)-Mu)
  Quadratic <- colSums(Eta*Eta)

  cur.alpha <- Alpha
  cur.beta <- Beta
  for(i in 2:max.iter){
    # -- update alpha
    routine.alpha <- optim(cur.alpha[i-1], method = "Brent", fn = function(a){
      -(
        a*(n*log(cur.beta[i-1]) - sum(log(Quadratic/2 + cur.beta[i-1])))-n*(lgamma(a)-lgamma(p/2+a))
      )},
      gr = function(a){
        -(
          (n*log(cur.beta[i-1]) - sum(log(Quadratic/2 + cur.beta[i-1])))-n*(digamma(a)-digamma(p/2+a))
        )
      },
      lower = 1+1e-4, upper = 10^6, control = list(maxit = 1000))
    if(routine.alpha$convergence != 0){
      stop("Convergence error in phi!")
    }
    cur.alpha[i] <- routine.alpha$par

    # --update beta
    routine.beta <- optim(cur.beta[i-1], method = "Brent", function(b){
      -(
        n*cur.alpha[i]*log(b)-(p/2+cur.alpha[i])*sum(log(Quadratic/2+b))
      )},
      gr = function(b){
        -(
          n*cur.alpha[i]/b-(p/2+cur.alpha[i])*sum(1/(Quadratic/2+b))
        )
      },
      lower = 1e-4, upper = 10^6, control = list(maxit = 1000))
    if(routine.beta$convergence != 0){
      stop("Convergence error in phi!")
    }
    cur.beta[i] <- routine.beta$par

    if((abs((cur.alpha[i]-cur.alpha[i-1])/cur.alpha[i]) < conv.criterion$tollerance) & (abs((cur.beta[i]-cur.beta[i-1])/cur.beta[i]) < conv.criterion$tollerance)){
      conv.counter <- conv.counter + 1
      if(conv.counter == conv.criterion$iterations){
        cat("Converged\n")
        break
      }
    }
    else {
      conv.counter <- 0
    }
  }

  return(list(alpha = cur.alpha[i], beta = cur.beta[i], Niter = i))
}


