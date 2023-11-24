#' Update Alpha and Beta estimation
#'
#' @param x an object of class `spartacus`.
#' @param n.neighbors the number of nearest neighbors for fitting the nearest-neighbor Gaussian process (NNGP) model. The default value is 20. ??Ha senso lasciare la libertà di impostarlo, forse no??
#' @param gene.names a vector of length `nrow(x$x)`, containing the gene names. If `NULL`, the function takes `gene.names = row.names(x$x)`.
#'
#' @return a `spartacus` object with the updated estimate of Alpha and Beta
#' @export
#'
updateAlphaBeta <- function(x, n.neighbors = 20, gene.names = NULL){
  if(class(x) != "spartacus") stop("x is not a spartacus object")
  X <- x$x
  coordinates <- x$coordinates
  alpha0 <- x$alpha
  tau0 <- x$tau
  Mu <- x$mu
  Delta <- x$delta
  Beta <- x$beta
  Phi <- x$phi
  Cs <- x$Cs
  Ds <- x$Ds
  if(is.null(gene.names)) gene.names <- row.names(X)
  if(is.null(gene.names)) gene.names <- 1:nrow(X)

  K <- ifelse(is.vector(Mu), 1, nrow(Mu))
  R <- length(Phi)

  table.Ds <- as.vector(table(Ds))
  if(length(unique(table.Ds)) == 1){
    iNN <- array(NA, dim = c(table.Ds[1], n.neighbors+1, length(table.Ds)))
    iNN <- lapply(1:dim(iNN)[3], function(x) iNN[,,x])
  } else{
    iNN <- mapply(function(n.row, n.col) matrix(NA, n.row, n.col), table.Ds, n.neighbors+1)
  }
  for(r in 1:R){
    ordering <- orderMaxMinFast(coordinates[Ds == r,], numpropose = nrow(coordinates[Ds == r,]))
    coordinates[Ds == r,] <- coordinates[Ds == r,][ordering,]
    X[, Ds == r] <- X[,Ds == r][,ordering]

    iNN[[r]] <- GpGp::find_ordered_nn(coordinates[Ds == r,], m = n.neighbors)
  }


  Alpha.new <- matrix(NA, K, R)
  Beta.new <- matrix(NA, K, R)
  for(k in 1:K){
    for(r in 1:R){
      Linv <- computeLinv(covparms = c(1, Phi[r], Delta[k,r]),
                          covfun_name = "exponential_isotropic",
                          locs.ordered = coordinates[Ds == r,],
                          iNN = iNN[[r]])
      Eta <- Linv%*%(t(X[Cs == k, Ds == r])-Mu[k,r])
      Qt <- colSums(Eta*Eta)

      # -- update Alpha and beta estimates
      estimation.AlphaBeta <- updateAlphaBeta.CoCluster(x = x[Cs == k, Ds == r],
                                                        locs.ordered = coordinates[Ds ==r, ],
                                                        Linv = Linv,
                                                        iNN = iNN[[r]],
                                                        Alpha = alpha0,
                                                        Tau = tau0,
                                                        Mu = Mu[k,r],
                                                        Delta = Delta[k,r],
                                                        Beta = Beta[k,r],
                                                        Phi = Phi[r],
                                                        max.iter = 1000)
      Alpha.new[k,r] <- estimation.AlphaBeta$alpha
      Beta.new[k,r] <- estimation.AlphaBeta$beta
    }
  }

  # ?? Sovrascrivo i vecchi valori o aggiungo nuovi parametri??
  x$Alpha.new = Alpha.new
  x$Beta.new = Beta.new
  return(x)
}
