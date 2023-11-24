#' Compute the gene-specific variances
#'
#' This function returns the summary statistics for the random variances estimated through SpaRTaCUS within the `K*R` blocks.
#'
#' @param x x an object of class `spartacus`.
#' @param x an object of class `spartaco`;
#' @param gene.names a vector of length `nrow(x$x)`, containing the gene names. If `NULL`, the function takes `gene.names = row.names(x$x)`.
#' @param MCMC.interval a list of length 2 containing the options for computing the high posterior density intervals of the gene specific variances. The  level is given by `MCMC.interval$level` (default `0.95`).
#' Intervals are computed via MCMC; the number of simulated draws must is passed to `MCMC.interval$level` (default is `10^4`).
#' If `NULL`, the intervals are computed as mean + c(-1.96,1.96)*sd, and the function returns a warning.
#'
#' @return An object of class `spartacus.genes` containing five elements:
#' - `Expectation`: a data frame of dimension `nrow(x$x)` x `R` containing the expected values of the random effects within the `R` spot clusters;
#' - `Variance`: a data frame of dimension `nrow(x$x)` x `R` containing the variance of the random effects within the `R` spot clusters;
#' - `HPD.left`: a data frame of dimension `nrow(x$x)` x `R` containing the left HPD interval of the random effects within the `R` spot clusters;
#' - `HPD.right` a data frame of dimension `nrow(x$x)` x `R` containing the right HPD interval of the random effects within the `R` spot clusters;
#' - `Cs`: the row clustering labels (the same contained into the object `x`).
#' - `Alpha.new`: the updated estiamte of the Alpha parameter of the SpaRTaCUS model.
#' - `Beta.new`: the updated estiamte of the Beta parameter of the SpaRTaCUS model.
#'
geneVariances.approx <- function(x, gene.names = NULL, MCMC.interval = list(level = .95, R = 10^4)){
  if(class(x) != "spartacus") stop("x is not a spartacus object")
  X <- x$x
  coordinates <- x$coord
  n.neighbors <- x$n.neighbors
  Alpha <- x$alpha
  Tau <- x$tau
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
  Alpha.post <- Beta.post <- Expected.post <- Variance.post <- data.frame(matrix(0, nrow(X), R))
  rownames(Alpha.post) <- rownames(Beta.post) <- rownames(Expected.post) <-
    rownames(Variance.post) <- gene.names
  colnames(Alpha.post) <- colnames(Beta.post) <- colnames(Expected.post) <-
    colnames(Variance.post) <- paste("r =",1:R)

  Inter.low <- Inter.up <- data.frame(matrix(0, nrow(X), R))
  rownames(Inter.low) <- rownames(Inter.up) <- gene.names
  colnames(Inter.low) <- colnames(Inter.up) <- paste("r =",1:R)


  for(k in 1:K){
    for(r in 1:R){
      Linv <- computeLinv(covparms = c(1, Phi[r], Delta[k,r]),
                          covfun_name = "exponential_isotropic",
                          locs.ordered = coordinates[Ds == r,],
                          iNN = iNN[[r]])
      Eta <- Linv%*%(t(X[Cs == k, Ds == r])-Mu[k,r])
      Qt <- colSums(Eta*Eta)


      # -- update Alpha and beta estimates
      estimation.AlphaBeta <- updateAlphaBeta.CoCluster(x = X[Cs == k, Ds == r],
                                                        locs.ordered = coordinates[Ds == r, ],
                                                        Linv = Linv,
                                                        Alpha = Alpha,
                                                        Tau = Tau,
                                                        Mu = Mu[k,r],
                                                        Delta = Delta[k,r],
                                                        Beta = Beta[k,r],
                                                        Phi = Phi[r],
                                                        max.iter = 10000,
                                                        conv.criterion = list(tollerance = 1e-3, iterations = 10))
      Alpha.new[k,r] <- estimation.AlphaBeta$alpha[length(estimation.AlphaBeta$alpha)]
      Beta.new[k,r] <- estimation.AlphaBeta$beta[length(estimation.AlphaBeta$beta)]

      Alpha.post[Cs == k,r] <- table.Ds[r]/2 + Alpha.new[k,r]
      Beta.post[Cs == k,r] <- Qt/2 + Beta.new[k,r]
      Expected.post[Cs == k,r] <- ifelse(Alpha.post[Cs == k,r]>1, Beta.post[Cs == k,r]/(Alpha.post[Cs == k,r]-1), Inf)
      Variance.post[Cs == k,r] <- ifelse(Alpha.post[Cs == k,r]>2, Beta.post[Cs == k,r]^2/((Alpha.post[Cs == k,r]-1)^2 * (Alpha.post[Cs == k,r]-2)), Inf)

      if(!is.null(MCMC.interval)){
        for(i in 1:sum(Cs == k)){
          gen <- rinvgamma(MCMC.interval$R, shape = Alpha.post[Cs == k,r][i], rate = Beta.post[Cs == k,r][i])
          interv <- HPDinterval(obj = mcmc(gen), prob = MCMC.interval$level)
          Inter.low[Cs == k,r][i] <- interv[1]
          Inter.up[Cs == k,r][i] <- interv[2]
        }
      } else {
        Inter.low[,r] <- Expected.post[,r] - 2*sqrt(Variance.post[,r])
        Inter.up[,r] <- Expected.post[,r] + 2*sqrt(Variance.post[,r])
      }
    }
  }
  if(is.null(MCMC.interval)) warning("the HPD intervals returned are approximated as mean + c(-1.96,1.96)*sd")
  output <- list(Expectation = Expected.post,
                 Variance = Variance.post,
                 HPD.left = Inter.low,
                 HPD.right = Inter.up,
                 Cs = Cs,
                 Alpha.new = Alpha.new,
                 Beta.new = Beta.new)
  class(output) <- "spartaco.genes"
  return(output)
}
