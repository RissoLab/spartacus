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
#' @param spot.labels a vector containing the labels of column cluster.
#' @param n.neighbors the number of nearest neighbors for fitting the nearest-neighbor Gaussian process (NNGP) model. The default value is 20. ??Da aggiungere commento sulla scelta di m??
#' @param K the number of row clusters (only when `input.values == NULL`);
#' @param Alpha the constraint for the Alpha parameters. ??Da aggiungere nei dettagli le caratteristiche del vincolo, se lasciamo la libertà di farlo??
#' @param Tau the constraint for the Tau parameters. ??Da aggiungere nei dettagli le caratteristiche del vincolo, se lasciamo la libertà di farlo??
#' @param compute.uncertainty if `TRUE` (default), it computes the clustering uncertainty of the rows and of the columns.
#' @param nstart the number of parallel runs of the estimation algorithm.
#' @param max.iter the maximum number of iterations the estimation algorithm is run.
#' @param estimate.iterations the maximum number of iterations within each M Step.
#' @param conv.criterion a list containing the parameters that define a converge criterion (see **Details**).
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
                               spot.labels,
                               approximated = TRUE,
                               n.neighbors = 20,
                               Alpha = 10,
                               Tau = 1,
                               compute.uncertainty = TRUE,
                               nstart = 5,
                               max.iter = 1000,
                               estimate.iterations = 100,
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
      spartacus.internal(data = x, coordinates = coordinates,
                         K = K, R = R, unsupervised = unsupervised, spot.labels = spot.labels,
                         n.neighbors = n.neighbors, Alpha = Alpha, Tau = Tau,
                         max.iter = max.iter, conv.criterion = conv.criterion,
                         verbose = F), future.seed = NULL
    )
  }
  else{
    # da inserire codice di spartaco con possibilità di specificare le etichette di spot
  }


  if(K == 1) compute.uncertainty = F
  output <- combineSpartacus(results, compute.uncertainty = compute.uncertainty)
  return(output)
}
