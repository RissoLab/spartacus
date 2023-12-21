main.approx <- function(x,
                        coordinates,
                        K,
                        R = NULL,
                        unsupervised,
                        spot.labels,
                        n.neighbors,
                        Alpha,
                        Tau,
                        input.values = NULL,
                        max.iter = 10^3,
                        Metropolis.iterations = 150,
                        conv.criterion = NULL,
                        save.options = NULL,
                        verbose = F
){
  if(verbose == T) cat("Parameters and labels allocation/ \n")
  Dist <- as.matrix(stats::dist(coordinates))
  R <- ifelse(unsupervised, R, length(unique(spot.labels)))
  if(is.null(input.values)){
    if(is.null(spot.labels)){
      cur.Ds <- best.Ds <- sample(1:R, size = ncol(x), replace = T)
    } else{
      cur.Ds <- as.numeric(spot.labels)
    }
    cur.Cs <- best.Cs <- sample(1:K, size = nrow(x), replace = T)
    cur.phi <- best.phi <- runif(R, 1, 5)
    cur.mu <- best.mu <- matrix(runif(K*R, 1, 10), K, R)
    cur.delta <- best.delta <- matrix(runif(K * R, 1e-7, 10), K, R)
    cur.beta <- best.beta <- matrix(runif(K*R, 1, 3), K, R)
  } else {
    if(is.null(spot.labels)){
      cur.Ds <- best.Ds <- input.values$Ds
    } else{
      cur.Ds <- as.numeric(spot.labels)
    }
    cur.Cs <- best.Cs <- input.values$Cs
    cur.mu <- best.mu <- input.values$mu
    cur.phi <- best.phi <- input.values$phi
    cur.delta <- best.delta <- input.values$delta
    cur.beta <- best.beta <- input.values$beta
  }
  goodK <- 1:K
  goodR <- 1:R

  if(verbose == T) cat("Inizialize iNN.full matrix/ \n")
  ordering <- orderMaxMinFast(coordinates, numpropose = nrow(coordinates))
  x <- x[,ordering]
  coordinates <- coordinates[ordering,]
  cur.Ds <- cur.Ds[ordering]
  iNN.full <- GpGp::find_ordered_nn(coordinates, m = nrow(coordinates))

  if(verbose == T) cat("Inizialize iNN matrix/ \n")
  iNN <- compute_iNN(Ds = cur.Ds, iNN.full = iNN.full, n.neighbors = n.neighbors)

  ll <- rep(-1e+40, max.iter)
  logL.values.approx <- matrix(0, K, R)
  i <- 1
  if(!is.null(conv.criterion)) counter.conv <- 0
  while(T){
    if(i == max.iter) break
    i <- i + 1

    #---M Step
    if(verbose == T) cat("iteration: ", i, " M Step/ ")
    for(r in goodR){
      for(k in goodK){
        estimation.parameters <- estimate.CoCluster.Parameters.marginal.approx(x = x[cur.Cs == k, cur.Ds == r],
                                                                               locs.ordered = coordinates[cur.Ds == r, ],
                                                                               iNN = iNN[[r]],
                                                                               Alpha = Alpha,
                                                                               Tau = Tau,
                                                                               Phi = cur.phi[r],
                                                                               mu0 = cur.mu[k,r],
                                                                               delta0 = cur.delta[k,r],
                                                                               beta0 = cur.beta[k,r])
        cur.mu[k,r] <- estimation.parameters$mu
        cur.delta[k,r] <- estimation.parameters$delta
        cur.beta[k,r] <- estimation.parameters$beta
      }
      cur.phi[r] <- estimatePhi.approx(x = x[,cur.Ds == r],
                                       locs.ordered = coordinates[cur.Ds == r, ],
                                       iNN = iNN[[r]],
                                       Cs = cur.Cs,
                                       Alpha = Alpha,
                                       Tau = Tau,
                                       Mu = cur.mu[,r],
                                       Delta = cur.delta[,r],
                                       Beta = cur.beta[,r],
                                       phi.old = cur.phi[r])
    }

    #---CE Step
    if(verbose == T) cat("CE Step/ ")
    cur.cs <- RowClustering.approx(x = x, Ds = cur.Ds, coordinates = coordinates, iNN = iNN,
                                   Tau = Tau, Alpha = Alpha, Mu = cur.mu, Delta = cur.delta, Beta = cur.beta, Phi = cur.phi)
    cur.Cs <- cur.cs$allocation

    #---SE Step
    if(verbose == T) cat("SE Step/\n")
    if(unsupervised){
      if(verbose == "full") cat("SE Step/")
      cur.ds <- MetropolisAllocation(x = x, Cs = cur.Cs, Ds = cur.Ds, Mu = cur.mu, Tau = Tau, Alpha = Alpha, Beta = cur.beta,
                                     Phi = cur.phi, Delta = cur.delta, coordinates = coordinates,
                                     iNN.full = iNN.full, iNN = iNN, n.neighbors = n.neighbors, maxit = Metropolis.iterations)
      cur.Ds <- cur.ds$Ds
      iNN <- compute_iNN(Ds = cur.Ds, iNN.full = iNN.full, n.neighbors = n.neighbors)
    }
    # se unsupervised questa parte può non essere fatta
    goodK <- sort(unique(cur.Cs))
    goodR <- sort(unique(cur.Ds))
    for(r in goodR){
      for(k in goodK){
        logL.values.approx[k,r] <- logL.Cocluster.approx(x = x[cur.Cs == k, cur.Ds == r],
                                                         locs.ordered = coordinates[cur.Ds == r,],
                                                         iNN = iNN[[r]],
                                                         Alpha = Alpha,
                                                         Tau = Tau,
                                                         Mu = cur.mu[k,r],
                                                         Beta = cur.beta[k,r],
                                                         Phi = cur.phi[r],
                                                         Delta = cur.delta[k,r])
      }
    }
    ll[i] <- sum(logL.values.approx)

    if(ll[i] == max(ll)){
      best.mu <- cur.mu
      best.delta <- cur.delta
      best.beta <- cur.beta
      best.phi <- cur.phi
      best.Cs <- cur.Cs
      best.Ds <- cur.Ds
    }

    if(!is.null(conv.criterion)){
      if(round(ll[i],6) >= round(ll[i-1],6) & (round(ll[i] - ll[i-1],6) < conv.criterion$epsilon)){
        counter.conv <- counter.conv + 1
        if(counter.conv == conv.criterion$iterations){
          cat("Converged\n")
          break
        }
      } else {
        counter.conv <- 0
      }
    }

    # Save the result in the given location after save.options$after iterations
    if(!is.null(save.options)){
      if(i %% save.options$after == 0){
        ICL <- max(ll) - nrow(x)*K - ncol(x)*R - .5*(4*K*R+R)*log(nrow(x) * ncol(x))
        results <- list(
          alpha = Alpha,
          tau = Tau,
          mu = best.mu,
          delta = best.delta,
          beta = best.beta,
          phi = best.phi,
          Cs = best.Cs,
          Ds = best.Ds[order(ordering)],
          logL = ll[c(2:i)],
          ICL = ICL,
          x = x[,order(ordering)],
          coordinates = coordinates[order(ordering),],
          n.neighbors = n.neighbors,
          unsupervised = unsuperviseda
        )
        class(results) <- "spartacus"
        save(results, file = save.options$file.name)
      }
    }
  }

  ICL <- max(ll) - nrow(x)*log(K) - ncol(x)*log(R) - .5*(4*K*R+R)*log(nrow(x) * ncol(x))

  # Save the result in the given location
  if(!is.null(save.options)){
    results <- list(
      alpha = Alpha,
      tau = Tau,
      mu = best.mu,
      delta = best.delta,
      beta = best.beta,
      phi = best.phi,
      Cs = best.Cs,
      Ds = best.Ds[order(ordering)],
      logL = ll[c(2:i)],
      ICL = ICL,
      x = x[,order(ordering)],
      coordinates = coordinates[order(ordering),],
      n.neighbors = n.neighbors,
      unsupervised = unsupervised
    )
    class(results) <- "spartacus"
    save(results, file = save.options$file.name)
  }
  results <- list(
    alpha = Alpha,
    tau = Tau,
    mu = best.mu,
    delta = best.delta,
    beta = best.beta,
    phi = best.phi,
    Cs = best.Cs,
    Ds = best.Ds[order(ordering)],
    logL = ll[c(2:i)],
    ICL = ICL,
    x = x[,order(ordering)],
    coordinates = coordinates[order(ordering),],
    n.neighbors = n.neighbors,
    unsupervised = unsupervised
  )
  class(results) <- "spartacus"
  return(results)
}

