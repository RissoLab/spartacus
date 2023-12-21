MetropolisAllocation <- function(x,
                                 Cs, Ds,
                                 Mu, Tau, Alpha, Beta, Phi, Delta,
                                 coordinates, iNN.full, iNN, n.neighbors,
                                 maxit = 150,
                                 min.obs = 3,
                                 prob.m = c(.7, .2, .1)){
  K <- ifelse(is.vector(Mu), 1, nrow(Mu))
  R <- length(Phi)
  if(is.vector(Mu)) Mu <- matrix(Mu, K, R)
  if(is.vector(Tau)) Tau <- matrix(Tau, K, R)
  if(is.vector(Delta)) Delta <- matrix(Delta, K, R)
  if(is.vector(Alpha)) Alpha <- matrix(Alpha, K, R)
  if(is.vector(Beta)) Beta <- matrix(Beta, K, R)
  D <- as.vector(table(Ds))
  accepted <- F
  accepted <- m.list <- numeric(maxit)
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
  logL.den <- sum(logL.values)
  for(iter in 1:maxit){
    # set the number of columns to reallocate at step iter
    m <- m.list[iter] <- sample(1:length(prob.m),1, prob = prob.m)
    move <- sample(c("M1","M2"), size = 1)
    if(move == "M1"){
      gr.start <- sample(1:R, 1) # select one clustering where to remove m columns
      gr.end <- sample(setdiff(1:R, gr.start), 1) # select one clustering where to add m columns
      j <- sample(1:D[gr.start], m, replace = F) # select the m columns to switch labels
      Ds.star <- Ds
      Ds.star[which(Ds == gr.start)][j] <- gr.end
      to.be.changed <- unique(c(gr.start, gr.end)) # update logL in those co-cluster
      logL.values.star <- logL.values
      iNN.star <- compute_iNN(Ds = Ds.star, iNN.full = iNN.full, n.neighbors = n.neighbors, which.cluster = to.be.changed)
      for(index in 1:length(to.be.changed)){
        r = to.be.changed[index]
        logL.values.star[r] <- 0
        for(k in goodK){
          logL.values.star[r] <- logL.values.star[r] + logL.Cocluster.approx(x = x[Cs == k, Ds.star == r],
                                                                             Mu = Mu[k,r],
                                                                             Tau = Tau[k,r],
                                                                             Alpha = Alpha[k,r],
                                                                             Beta = Beta[k,r],
                                                                             Phi = Phi[r],
                                                                             Delta = Delta[k,r],
                                                                             locs.ordered = coordinates[Ds.star == r, ],
                                                                             iNN = iNN.star[[index]])
        }
      }
      log.proposal.num <- sum(log(D[gr.start]-0:(m-1)))
      log.proposal.den <- sum(log(D[gr.end]+1:m))
      logL.num <- sum(logL.values.star)
      A <- exp(logL.num + log.proposal.num - logL.den - log.proposal.den)*all(as.vector(table(Ds.star)) >= min.obs)
      if(runif(1) <= A){
        Ds <- Ds.star
        D <- as.vector(table(Ds))
        logL.values <- logL.values.star
        logL.den <- logL.num
        accepted[iter] <- 1
        #iNN[to.be.changed] <- iNN.star
      }
    }
    if(move == "M2"){
      gr.start <- sample(1:R, m, replace = T) # seleziono i cluster da cui rimuovere gli spot
      gr.end <- sapply(1:m, function(k) sample(setdiff(1:R, gr.start[k]), 1)) # seleziono i cluster in cui inserire gli spot
      q1r <- sapply(1:R, function(r) sum(gr.start == r)) # vettore con numero di spot da rimuovere per ogni cluster di spot
      q2r <- sapply(1:R, function(r) sum(gr.end == r)) # vettore con numero di spot da inserire per ogni cluster di spot
      j <- sapply(1:R, function(r) ifelse(q1r[r] != 0, return(sample(1:D[r], q1r[r], replace = F)), 0)) # seleziono gli spot da rimuovere per ogni cluster di spot
      Ds.star <- Ds
      for(r in which(q1r != 0)) Ds.star[which(Ds == r)][j[[r]]] <- gr.end[gr.start == r] # assegno nuova etichetta
      to.be.changed <- unique(c(gr.start, gr.end))
      logL.values.star <- logL.values
      iNN.star <- compute_iNN(Ds = Ds.star, iNN.full = iNN.full, n.neighbors = n.neighbors, which.cluster = to.be.changed)
      for(index in 1:length(to.be.changed)){
        r = to.be.changed[index]
        logL.values.star[r] <- 0
        for(k in goodK){
          logL.values.star[r] <- logL.values.star[r] + logL.Cocluster.approx(x = x[Cs == k,Ds.star == r],
                                                                             Mu = Mu[k,r],
                                                                             Tau = Tau[k,r],
                                                                             Alpha = Alpha[k,r],
                                                                             Beta = Beta[k,r],
                                                                             Phi = Phi[r],
                                                                             Delta = Delta[k,r],
                                                                             locs.ordered = coordinates[Ds.star == r, ],
                                                                             iNN = iNN.star[[index]])
        }
      }
      log.proposal.num <- sum(log(factorial(q2r)))+sum(sapply(which(q1r != 0), function(r) sum(log(D[r]+0:(q1r[r]-1)))))
      log.proposal.den <- sum(log(factorial(q1r)))+sum(sapply(which(q2r != 0), function(r) sum(log(D[r]-q1r[r]+1:q2r[r]))))
      logL.num <- sum(logL.values.star)
      A <- exp(logL.num + log.proposal.num - logL.den - log.proposal.den)*all(as.vector(table(Ds.star)) >= min.obs)
      if(runif(1) <= A){
        Ds <- Ds.star
        D <- as.vector(table(Ds))
        logL.values <- logL.values.star
        logL.den <- logL.num
        accepted[iter] <- 1
        #iNN[to.be.changed] <- iNN.star
      }
    }
  }
  results <- list(Ds = Ds, logL = sum(logL.values), logL.values = logL.values)
  return(results)
}
