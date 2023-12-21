#' SpaRTaCUS with Column Cluster assigned
#'
#' This function returns the estimated model parameters and the row-clustering labels.
#'
#' @import SpatialExperiment
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
#' @param max.iter the maximum number of iterations the estimation algorithm is run.
#' @param conv.criterion a list containing the parameters that define a converge criterion (see **Details**).
#' @param save.options a list for specifying the saving parameters (see **Details**).
#' @param seed set the internal seed of the function.
#' @param input.values starting value for parameter estimation
#' @param verbose progress of estimation
#' @param R numeber of spot cluster
#' @param unsupervised if true unspervised clustering is performed
#'
#' @return An object of class `spartacus` with the parameter estimates, the clustering labels, the log-likelihood value at each iteration and the data, the ICL, the data matrix and the coordinates matrix.
#'
#' @details
#'
#' The algorithm can be initiated from a given set of starting values. To do so, `input.values` receives a list of the form
#' `list(mu, tau, xi, alpha, beta, phi, Cs)`, where:
#' - `mu`, `delta` and `beta` are `K` x `R` matrices;
#' - `phi` is a vector of length `R`;
#' - `Cs` is a vector of length `nrow(data)` containing the row clustering labels;
#'
#' If the algorithm is initiated from some starting values,  `K` is set automatically according to the input values.
#' If an object of class `spartacus` is passed to `input.values`, the estimation starts from the final estimate of the previous run (see **Examples**).
#'
#' If `conv.criterion == NULL`, the algorithm is stopped after `max.iter` itereations, otherwise it is stopped when the increment of the log-likelihood is smaller than a certain threshold `conv.criterion$epsilon` for `conv.criterion$iterations` times in a row.
#'
#' The function allows also to save the results even if the estimation is not completed. To do so, `save.options` receives a list of two parameters:
#' `after` gives the number of iterations after which the results are saved, `file.name` contains the path where the results are saved.
#'
#' If `verbose == "full"`, the on-going estimation procedure is displayed. If `verbose == "progress"`, a dynamic progress bar will display the percentage of iterations completed.
#' If `verbose == F`, then nothing is displayed in console.

spartacus.internal <- function(data,
                              assay = NULL,
                              coordinates = NULL,
                              K = K, R = R,
                              unsupervised,
                              spot.labels,
                              n.neighbors = 20,
                              Alpha = 10,
                              Tau = 1,
                              input.values = NULL,
                              max.iter = 1000,
                              conv.criterion = list(iterations = 10, epsilon = 1e-4),
                              Metropolis.iterations = 150,
                              verbose = F,
                              save.options = NULL,
                              seed = NULL
) {

  if(class(data)[1] == "SpatialExperiment"){
    if(is.numeric(assay)) which.assay <- assay
    else which.assay <- which(names(data@assays@data) == assay)
    x <- as.matrix(data@assays@data[[which.assay]])
    row.names(x) <- rowData(data)$gene_name
    coordinates <- as.matrix(spatialCoords(data))
  } else {
    x <- data
  }

  set.seed(seed = seed)

  if(!is.null(input.values)){
    K <- nrow(input.values$mu)
  }

  main.approx(x = x, coordinates = coordinates,
              K = K, R = R, unsupervised = unsupervised, spot.labels = spot.labels,
              n.neighbors = n.neighbors, Alpha = Alpha, Tau = Tau,
              input.values = input.values,
              max.iter = max.iter, Metropolis.iterations = Metropolis.iterations, conv.criterion = conv.criterion,
              verbose = verbose, save.options = save.options
  )
}


