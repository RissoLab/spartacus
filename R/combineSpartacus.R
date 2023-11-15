#' Manually combine multiple runs of SpaRTaCUS
#'
#' This function combines multiple runs of SpaRTaCUS with the same values of K and R.

#' @import ggplot2
#' @export
#'
#' @param x either a list of `spartacus` objects or a vector containing the names of the files to be combined;
#' @param search.dir the directory path where the file names given by `x` are searched. If `NULL`, the file paths must be specified as part of the file names is `x`;
#' @param compute.uncertainty if `TRUE` (default), it computes the clustering uncertainty of the rows.
#'
#' @return an object of class `spartacus` with the parameter estimates and the clustering labels given by the best fit across the ones in `x`. If `compute.uncertainty == T`, the clustering uncertainty measures are returned.
#'
#' @details If each element of `x` is an object of class `spartacus`, then the function combines them into a unique spartaco object.
#' If the `spartaco` objects are stored into different files, the vector `x` receives the names of the files,
#' and `search.dir` is the path of the directory where the files are stored.

combineSpartacus <- function(x = NULL, search.dir = NULL, compute.uncertainty = T){
  loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
  }
  if(is.list(x)) results <- x else {
    results <- list()
    for(i in 1:length(x)){
      file.name <- x[i]
      if(!is.null(search.dir)){
        if(substr(search.dir, nchar(search.dir),nchar(search.dir)) == "/")
          file.name <- paste(search.dir, file.name, sep="") else
            file.name <- paste(search.dir, file.name, sep="/")
      }
      results[[i]] <- loadRData(file.name)
    }
  }
  likelihoods <- list()
  maxi <- rep(-Inf, length(results))

  final <- list()
  for(i in 1:length(results)){
    likelihoods[[i]] <- results[[i]]$logL
    maxi[i] <- max(likelihoods[[i]])
    if(which.max(maxi) == i){
      final$mu <- results[[i]]$mu
      final$tau <- results[[i]]$tau
      final$delta <- results[[i]]$delta
      final$stn.ratio <- results[[i]]$stn.ratio
      final$alpha <- results[[i]]$alpha
      final$beta <- results[[i]]$beta
      final$phi <- results[[i]]$phi
      final$Cs <- results[[i]]$Cs
      final$Ds <- results[[i]]$Ds
      final$ICL <- results[[i]]$ICL
      final$logL <- results[[i]]$logL
      final$n.neighbors <- results[[i]]$n.neighbors
    }
  }
  final$max.logL <- maxi
  final$x <- results[[i]]$x
  final$coordinates <- results[[i]]$coordinates

  if(compute.uncertainty){
    if(nrow(final$mu) == 1) return(-1)
    best.j <- which.max(final$max.logL)
    j.to.invest <- setdiff(1:length(final$max.logL), best.j)
    final$cluster.discr <- list()
    # evaluate the discrepancy across the rows
    cat("Computing the row discrepancy...\n")
    CERs <- matrix(0, nrow(final$mu), length(j.to.invest))
    sapply(1:nrow(final$mu), function(k){
      reference <- as.numeric(final$Cs == k)
      sapply(1:length(j.to.invest), function(j){
        classif <- table(final$Cs, results[[j.to.invest[j]]]$Cs)[k,]
        k.j <- as.numeric(attr(classif, "names"))[which.max(classif)]
        comparison <- as.numeric(results[[j.to.invest[j]]]$Cs == k.j)
        CERs[k,j] <<- CER(reference = reference, estimate = comparison)
      })
    })
    w <- 1/(final$max.logL[best.j]-final$max.logL[-best.j])
    final$cluster.discr$rows <- as.vector((CERs %*% w)/sum(w))
  }

  # execution time
  final$exec.time <- numeric(length(x))
  sapply(1:length(x), function(i)
    if(is.null(results[[i]]$exec.time)) final$exec.time[i] <<- NA else
      final$exec.time[i] <<- results[[i]]$exec.time)

  ChainPlot <- data.frame(
    ranges = as.vector(unlist(likelihoods)),
    iterations = as.vector(unlist(sapply(1:length(x), function(j) 1:length(likelihoods[[j]])))),
    chains = as.factor(as.vector(unlist(sapply(1:length(x), function(j) rep(j, length(likelihoods[[j]]))))))
  )
  final$ll.plot <- ggplot(ChainPlot, aes(iterations, ranges, col = chains))+geom_line()+
    theme_classic()+theme(legend.position = "bottom")+
    labs(x = "iteration", y = "classification log-likelihood", color = "run")
  plot(final$ll.plot)
  # plot(1,1,
  #     xlim = c(1, max(sapply(1:length(likelihoods), function(l) length(likelihoods[[l]])))),
  #     ylim = c(min(ranges), max(ranges)),
  #     xlab = "Iteration", ylab = expression(logL))
  #for(i in 1:length(likelihoods))
  #    points(1:length(likelihoods[[i]]), likelihoods[[i]], col = i, pch = i)

  class(final) <- "spartacus"
  return(final)
}
