#' Multiple runs of SpaRTaCUS
#'
#' This function returns the estimated model parameters and the co-clustering labels obtained after running SpaRTaCUS multiple times (parallel options available).
#'
#' @import SpatialExperiment
#' @import future.apply
#' @import future
#' @import fields
#' @import Matrix
#' @import matrixStats
#' @import GpGp
#'
#' @export
#'
#' @param data either a `SpatialExperiment` object or a matrix containing the experiment.
#' @param assay if `class(data) == "SpatialExperiment"`, it takes either the name or the index of the assay to be used.
#' @param coordinates if `is.matrix(data)`, it takes the matrix of spatial coordinates of dimension `ncol(data)` x 2.
#' @param unsupervised set if the columns clustering has to be performed or not, the default is `TRUE`. If `FALSE`, spot.labels must be provided.
#' @param spot.labels a vector containing the labels of column cluster, the default is `NULL`.
#' @param approximated if `TRUE` the spartacus method is performed (default), else if `FALSE` the spartaco method is performed.
#' @param n.neighbors the number of nearest neighbors for fitting the nearest-neighbor Gaussian process (NNGP) model. The default value is 20. ??Da aggiungere commento sulla scelta di m??
#' @param K the number of row clusters (only when `input.values == NULL`);
#' @param R the number of column clusters (only when `input.values == NULL`);
#' @param Alpha the constraint for the Alpha parameters. ??Da aggiungere nei dettagli le caratteristiche del vincolo, se lasciamo la libertà di farlo??
#' @param Tau the constraint for the Tau parameters. ??Da aggiungere nei dettagli le caratteristiche del vincolo, se lasciamo la libertà di farlo??
#' @param nstart the number of parallel runs of the estimation algorithm.
#' @param max.iter the maximum number of iterations the estimation algorithm is run.
#' @param conv.criterion a list containing the parameters that define a converge criterion (see **Details**).
#' @param Metropolis.iterations the number of iteration done in the column clustering estimation.
#'
#' @return An object of class `spartacus` with the parameter estimates, the row clustering labels, the log-likelihood value at each iteration, the ICL, the data matrix and the coordinates matrix, and the clustering uncertainty.
#'
#' @details This function allows to run the `spartacus` model starting from multiple starting points simultaneously and select the one with the highest likelihood.
#' It is possible to run this function using multiple cores; to do so, use the `multicore` function package `future` (see **Examples**).
#'
#' ??Da aggiornare??
#' @examples
#' library(spartacus)
#'
#' # First, create the data matrix:
#' n <- p <- 300
#' K <- R <- 3
#' x <- matrix(runif(n*p), n, p)
#' coordinates <- matrix(runif(2*p), p, 2)
#'
#' # Set the number of cores to be used for the computations. In this example, we use 3 cores.
#' future::plan(future::multisession(workers = 3))
#' output <- spartacus_multirun(data = x, coordinates = coordinates, K = K, R = R, max.iter = 1000)
#'

spartacus_multirun <- function(data,
                               assay = NULL,
                               coordinates = NULL,
                               K,
                               R = NULL,
                               unsupervised = TRUE,
                               spot.labels = NULL,
                               approximated = TRUE,
                               n.neighbors = 20,
                               Alpha = 10,
                               Tau = 1,
                               nstart = 5,
                               max.iter = 1000,
                               Metropolis.iterations = 150,
                               conv.criterion = list(iterations = 10, epsilon = 1e-4)
)
{
  if(class(data)[1] == "SpatialExperiment"){ # è il modo giusto di gestire il caso in cui class() è un vettore? Ad esempio se data è la matrice di espressione dei geni con i nomi di riga e colonna class è c("matrix", "array")
    if(is.numeric(assay)) which.assay <- assay
    else which.assay <- which(names(data@assays@data) == assay)
    x <- as.matrix(data@assays@data[[which.assay]])
    row.names(x) <- rowData(data)$gene_name
    coordinates <- as.matrix(spatialCoords(data))
  } else {
    x <- data
  }

  if(approximated){
    results <- future_lapply(1:nstart, FUN = function(l)
      spartacus.internal(data = x, assay = NULL, coordinates = coordinates, R = R, K = K, unsupervised = unsupervised, spot.labels = spot.labels,
                         max.iter = max.iter, Metropolis.iterations = Metropolis.iterations, Alpha = Alpha, Tau = Tau,
                         conv.criterion = conv.criterion, verbose = F, save.options = NULL, n.neighbors = n.neighbors)
      , future.seed = TRUE)
  }
  else{
    # da inserire codice di spartaco con possibilità di specificare le etichette di spot
  }


  compute.uncertainty.row <- ifelse(K == 1, F, T)
  compute.uncertainty.col <- ifelse(R == 1, F, T)
  output <- combineSpartacus(results, compute.uncertainty.row = compute.uncertainty.row, compute.uncertainty.col = compute.uncertainty.col)
  return(output)
}
