computeLinv <- function(covparms, covfun_name = "exponential_isotropic", locs.ordered, iNN, m){
  Linv <- GpGp::vecchia_Linv(covparms = covparms,
                       covfun_name = covfun_name,
                       locs = locs.ordered,
                       NNarray = iNN)

  # The function "vecchia_Linv" returns ONLY the entries of the inverse Cholesky factor of the covariance matrix implied by Vecchia's approximation
  # Here the full Linv matrix is made
  howMany = c(seq_len(m),rep(m+1, nrow(locs.ordered)-m)) # number of non-zero elements by row
  iRow = rep(1:nrow(locs.ordered), times = howMany)
  iCol = na.omit((as.vector(t(iNN)))) # number of non-zero elements by column
  L.elms = as.vector(t(Linv))[!is.na((as.vector(t(iNN))))]
  Linv <- sparseMatrix(i = iRow, j = iCol, x = L.elms, triangular = TRUE)

  return(Linv)
}
