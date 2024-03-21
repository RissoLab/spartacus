GibbsAllocation <- function(x,
                            Cs, Ds,
                            Mu, Tau, Alpha, Beta, Phi, Delta,
                            coordinates, iNN.full, iNN, n.neighbors,
                            maxit = 10,
                            min.obs = 3,
                            prob.m = c(.7, .2, .1)){
  K <- ifelse(is.vector(Mu), 1, nrow(Mu))
  R <- length(Phi)
  if(is.vector(Mu)) Mu <- matrix(Mu, K, R)
  if(is.vector(Tau)) Tau <- matrix(Tau, K, R)
  if(is.vector(Delta)) Delta <- matrix(Delta, K, R)
  if(is.vector(Alpha)) Alpha <- matrix(Alpha, K, R)
  if(is.vector(Beta)) Beta <- matrix(Beta, K, R)
  goodK <- sort(unique(Cs))
  goodR <- sort(unique(Ds))
  logL.values <- numeric(R)
  for(r in goodR){
    for(k in goodK){
      logL.values[r] <- logL.values[r] + logL.Cocluster.approx(x = x[Cs == k, Ds == r],
                                                               Mu = Mu[k,r],
                                                               Tau = Tau[k,r],
                                                               Alpha = Alpha[k,r],
                                                               Beta = Beta[k,r],
                                                               Phi = Phi[r],
                                                               Delta = Delta[k,r],
                                                               locs.ordered = coordinates[Ds == r,],
                                                               iNN = iNN[[r]])
    }
  }
  for(j in 1:ncol(x)){
    Ds.star <- Ds
    probabilities <- numeric(R)
    # inizializzo ll con i valori nella configurazione precedenete poi aggiorno solo per i co-cluster che cambiano
    ll <- matrix(logL.values, R, R, byrow = T)

    # move a single observation
    r.start <- Ds[j]
    Ds.star[j] <- R+1

    # starting cluster - update ll
    iNN.star <- compute_iNN(Ds = Ds.star, iNN.full = iNN.full, n.neighbors = n.neighbors, which.cluster = r.start)
    ll[-r.start,r.start] <- 0
    for(k in goodK){
      ll[-r.start,r.start] <- ll[-r.start,r.start] + logL.Cocluster.approx(x = x[Cs == k, Ds.star == r.start],
                                                                           Mu = Mu[k,r.start],
                                                                           Tau = Tau[k,r.start],
                                                                           Alpha = Alpha[k,r.start],
                                                                           Beta = Beta[k,r.start],
                                                                           Phi = Phi[r.start],
                                                                           Delta = Delta[k,r.start],
                                                                           locs.ordered = coordinates[Ds.star == r.start, ],
                                                                           iNN = iNN.star[[1]])
    }
    for(r.end in setdiff(1:R, r.start)){
      # receiving cluster
      Ds.star[j] <- r.end
      iNN.star <- compute_iNN(Ds = Ds.star, iNN.full = iNN.full, n.neighbors = n.neighbors, which.cluster = r.end)
      ll[r.end,r.end] <- 0
      for(k in goodK){
        ll[r.end,r.end] <- ll[r.end,r.end] + logL.Cocluster.approx(x = x[Cs == k, Ds.star == r.end],
                                                                   Mu = Mu[k,r.end],
                                                                   Tau = Tau[k,r.end],
                                                                   Alpha = Alpha[k,r.end],
                                                                   Beta = Beta[k,r.end],
                                                                   Phi = Phi[r.end],
                                                                   Delta = Delta[k,r.end],
                                                                   locs.ordered = coordinates[Ds.star == r.end, ],
                                                                   iNN = iNN.star[[1]])
      }
    }
    probabilities <- rowSums(ll)
    #--- stochastic allocation
    pp <- exp(probabilities-max(probabilities))
    pp <- pp/sum(pp)
    Ds[j] <- sample(1:R, prob = pp, size = 1)
    #Ds[j] <- which.max(probabilities)
    logL.values <- ll[Ds[j],]

    if(any(as.vector(table(Ds)) < 3)){
      results <- list(Ds = Ds, logL = sum(logL.values), logL.values = logL.values)
      return(results)
    }
  }

  results <- list(Ds = Ds, logL = sum(logL.values), logL.values = logL.values)
  return(results)
}



