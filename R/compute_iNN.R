compute_iNN <- function(Ds, iNN.full, n.neighbors, which.cluster = NULL){
  if(is.null(which.cluster)){
    which.cluster <- 1:length(unique(Ds))
    goodR <- 1:length(which.cluster)
    table.Ds <- as.vector(table(Ds))
  } else {
    table.Ds <- as.vector(table(Ds[Ds %in% which.cluster]))
    goodR <- 1:length(which.cluster)
  }
  n.neighbors.per.cluster <- sapply(table.Ds, function(s) min(s, n.neighbors-1))
  # capire se è meglio allocarla nel main
  iNN <- mapply(function(n.row, n.col) matrix(NA, n.row, n.col), table.Ds, n.neighbors.per.cluster+1, SIMPLIFY = F)
  for(r in goodR){
    iNN.full2 <- iNN.full[Ds == which.cluster[r],]
    iNN.full2[!(iNN.full2 %in% which(Ds == which.cluster[r]))] <- NA
    iNN[[r]] <- t(sapply(1:nrow(iNN.full2), function(i){
      temp <- iNN.full2[i,!is.na(iNN.full2[i,])]
      if(length(temp) >= n.neighbors.per.cluster[r]) temp <- temp[1:n.neighbors.per.cluster[r]] else
        temp <- c(temp, rep(NA,n.neighbors.per.cluster[r]-length(temp)))
      return(temp)
    }))
    coefs <- as.vector(iNN[[r]])
    ranks <- rank(unique(coefs), na.last = NA)
    iNN[[r]] <- matrix(ranks[match(coefs, unique(coefs))], nrow(iNN[[r]]), ncol(iNN[[r]]))
  }
  return(iNN)
}
