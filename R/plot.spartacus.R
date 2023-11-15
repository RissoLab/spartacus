#' Plot SpaRTaCo
#' This function returns the ggplots of the mean and the spatial signal-to-noise ratios from a SpaRTaCUS model.
#'
#' @import ggplot2
#' @export
#'
#' @param x a `spartacus` object;
#' @param type the type of plot you want to show.
#' - `1` displays the co-cluster mean levels;
#' - `2` displays the co-cluster spatial signal-to-noise ratios;
#' - `3` displays the map of the spots colored with respect to the estimated column clusters;
#' - `4` if `!is.null(k)`, it displays the average gene expression in the clusters `k` and `r`, otherwise it displays the expression of the gene given in `gene.name`.
#' @param gene.name  (used only when `type == 4`); it receives the name of the gene to display. If `!is.null(k)`, it displays the average expression of the gene clusters given by `k`.
#' @param k (used only when `type == 4`) the gene clusters to plot.
#' @param r (used only when `type == 4`) the spot clusters to plot.
#' @param manual.palette a vector of colors used when `type == 3`.
#' @param display.all.spots if `TRUE` (default) and `type == 4`, it displays the entire grid of spots.
#' @param range.mu set min and max values of mean's plot scale.
#' @param range.sRatio set min and max values of spatian signal to noise Ratio's plot scale.
#' @param ...
#'
#' @return The requested plot is displayed. In addition, if assigned to an object, it will return the `ggplot` object.
#'
plot.spartacus <- function(x, type = 1, gene.name = readline("gene name: "), k = NULL, r = 1:ncol(x$mu), manual.palette = NULL,
                          display.all.spots = T, range.mu = NULL, range.sRatio = NULL, ...){
  if(class(x) != "spartacus") stop("the input file is not a spartacus object")
  if(length(type) > 1) type <- 1
  K <- nrow(x$mu)
  R <- ncol(x$mu)
  k.lab <- 1:K
  r.lab <- 1:R
  gr <- expand.grid(k.lab,r.lab)
  gr$Mu <- as.vector(x$mu)
  gr$Ratio <- as.vector(1/x$delta)
  gr <- cbind(gr, expand.grid(as.vector(table(x$Cs))/nrow(x$x),as.vector(table(x$Ds))/ncol(x$x)))
  names(gr)[-c(3,4)] <- c("Y","X","height","width")

  prop.x <- as.vector(table(x$Ds))/ncol(x$x)
  prop.y <- as.vector(table(x$Cs))/nrow(x$x)
  xlim.left <- c(0, cumsum(prop.x)[-R])
  xlim.right <- cumsum(prop.x)
  ylim.left <- c(0, cumsum(prop.y)[-K])
  ylim.right <- cumsum(prop.y)

  if(is.null(range.mu)){
    mu.min <- min(gr$Mu)
    mu.max <- max(gr$Mu)
  }else{
    mu.min <- min(range.mu)
    mu.max <- max(range.mu)
  }

  if(is.null(range.sRatio)){
    sRatio.min <- min(gr$Ratio)
    sRatio.max <- max(gr$Ratio)
  }else{
    sRatio.min <- min(range.sRatio)
    sRatio.max <- max(range.sRatio)
  }
  # ---plot mu
  if(type == 1){
    Plots <- ggplot(gr, aes(xmin = as.vector(sapply(1:R, function(i) rep(xlim.left[i],K))),
                            xmax = as.vector(sapply(1:R, function(i) rep(xlim.right[i],K))),
                            ymin = rep(ylim.left,R),
                            ymax = rep(ylim.right,R),
                            fill = Mu)
    )+geom_rect()+theme_bw()+
      scale_fill_distiller(palette = "RdPu", limits = c(mu.min, mu.max),
                           labels = round(seq(mu.min, mu.max, length = 5),2),
                           breaks = seq(mu.min, mu.max, length = 5))+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text=element_text(size=18),
            axis.title=element_text(size=18),
            legend.text = element_text(size = 14),
            legend.title = element_text(size = 22),
            #legend.position = "bottom",
            plot.title = element_text(hjust = 0.5),
            title = element_text(size=18),
            plot.margin=grid::unit(c(3,2,3,2), "mm"))+
      labs(fill=expression(hat(mu)[kr]))+
      scale_x_continuous(breaks=(xlim.left+xlim.right)/2,
                         labels=paste("r =",1:R)
      )+
      scale_y_continuous(breaks=(ylim.left+ylim.right)/2,
                         labels=paste("k =",1:K)
      )


  }

  # ---plot tau/xi
  if(type == 2){
    Plots <- ggplot(gr, aes(xmin = as.vector(sapply(1:R, function(i) rep(xlim.left[i],K))),
                            xmax = as.vector(sapply(1:R, function(i) rep(xlim.right[i],K))),
                            ymin = rep(ylim.left,R),
                            ymax = rep(ylim.right,R),
                            fill = Ratio)
    )+geom_rect()+theme_bw()+
      viridis::scale_fill_viridis(discrete=FALSE, limits = c(sRatio.min, sRatio.max),
                                  labels = round(seq(sRatio.min, sRatio.max, length = 5),2),
                                  breaks = seq(sRatio.min, sRatio.max, length = 5))+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text=element_text(size=18),
            axis.title=element_text(size=18),
            legend.text = element_text(size = 14),
            legend.title = element_text(size = 22),
            #legend.position = "bottom",
            plot.title = element_text(hjust = 0.5),
            title = element_text(size=18),
            plot.margin=grid::unit(c(3,2,3,2), "mm"))+
      labs(fill=expression(hat(tau)[kr]/hat(xi)[kr]))+
      scale_x_continuous(breaks=(xlim.left+xlim.right)/2,
                         labels=paste("r =",1:R)
      )+
      scale_y_continuous(breaks=(ylim.left+ylim.right)/2,
                         labels=paste("k =",1:K)
      )


  }

  # ---plot spot clusters
  if(type == 3){
    # ---plot column clusters
        Coord <- data.frame(x = x$coordinates[,2], y = -x$coordinates[,1], z = as.factor(x$Ds))
        Coord$group <- as.factor(as.numeric(x$Ds %in% c(1,7:9)))
        Plots <- ggplot(Coord, aes(x, y, color = z))+
          geom_point(size = 3)+theme_bw()+
          labs(col = "")+
          labs(col = expression(D[r]))+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.x = element_blank(),#element_text(size=18),
                axis.text.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title=element_blank(),#element_text(size=18),
                legend.text = element_text(size = 18),
                legend.title = element_text(size = 22),
                plot.title = element_text(hjust = 0.5),
                #legend.position = "bottom",
                #legend.spacing.x = unit(0.3, 'cm'),
                title = element_text(size=18),
                plot.margin=grid::unit(c(3,2,3,2), "mm"))+
          geom_point(shape = 1,size = 3,colour = "black")
        #if(use.greys)
        #    Plots <- Plots + scale_color_grey(start = 0, end = .9) else
        Plots <- Plots + scale_fill_brewer( palette = "Set1")
  }

  # ---plot expressions
  if(type == 4){
    if(is.null(k)){
      if(gene.name != "all" & !(gene.name %in% row.names(x$x)))
        stop(paste("Gene",gene.name,"not found\n"))
      if(gene.name == "all") k <- 1:K}
    # ---plot sample means
    Coord <- data.frame(x = x$coordinates[,2], y = -x$coordinates[,1], z = as.factor(x$Ds))
    if(is.null(r)) r <- 1:R
    if(!display.all.spots){
      if(length(k) <= 1){
        if(length(k) == 0){
          x.bar <- x$x[which(row.names(x$x) == gene.name), which(x$Ds %in% r)]
        }
        if(length(k) == 1){
          x.bar <- colMeans(x$x[which(x$Cs == k), which(x$Ds %in% r)])
        }
        Plots <- ggplot(Coord[which(x$Ds %in% r),], aes(x, y, color = x.bar))+
          geom_point(size = 3)+theme_bw()+
          scale_fill_distiller(type = "seq", palette = "Spectral", direction = -1)+
          labs(col = "")+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.x = element_blank(),#element_text(size=18),
                axis.text.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title=element_blank(),#element_text(size=18),
                legend.text = element_text(size = 18),
                legend.title = element_text(size = 22),
                plot.title = element_text(hjust = 0.5),
                #legend.position = "bottom",
                #legend.spacing.x = unit(0.3, 'cm'),
                title = element_text(size=18),
                plot.margin=grid::unit(c(3,2,3,2), "mm"))+
          geom_point(shape = 1,size = 3,colour = "black")+
          ggtitle(label = ifelse(length(k) == 0, gene.name, paste("k =",k)))
      } else {
        Plots <- list()
        for(k.ind in 1:length(k)){
          Plots[[k.ind]] <- local({
            x.bar <- colMeans(x$x[x$Cs == k[k.ind], which(x$Ds %in% r)])
            ggplot(Coord[which(x$Ds %in% r),], aes(x, y, color = x.bar))+
              geom_point(size = 3)+theme_bw()+
              scale_color_distiller(type = "seq", palette = "Spectral", direction = -1)+
              labs(col = "")+
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.text.x = element_blank(),#element_text(size=18),
                    axis.text.y = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.title=element_blank(),#element_text(size=18),
                    legend.text = element_text(size = 18),
                    legend.title = element_text(size = 22),
                    plot.title = element_text(hjust = 0.5),
                    title = element_text(size=18),
                    plot.margin=grid::unit(c(3,2,3,2), "mm"))+
              geom_point(shape = 1,size = 3,colour = "black")+
              ggtitle(label = paste("k =",k[k.ind]))})
        }
      }
    } else {  # if display.all.spots == T
      if(length(k) <= 1){
        if(length(k) == 0){
          x.bar <- x$x[which(row.names(x$x) == gene.name), which(x$Ds %in% r)]
        }
        if(length(k) == 1){
          x.bar <- colMeans(x$x[which(x$Cs == k), which(x$Ds %in% r)])
        }
        Plots <- ggplot(Coord[-which(x$Ds %in% r),], aes(x, y))+
          theme_bw()+
          geom_point(data = Coord[-which(x$Ds %in% r),], mapping = aes(x, y), size = 3, fill = "white", colour = "gray74", shape = 21)+
          geom_point(data = Coord[which(x$Ds %in% r),], mapping = aes(x, y, fill = x.bar), colour = "white", size = 4, shape = 21)+
          scale_fill_distiller(type = "seq", palette = "Spectral", direction = -1)+
          labs(fill = "")+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.x = element_blank(),#element_text(size=18),
                axis.text.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title=element_blank(),#element_text(size=18),
                legend.text = element_text(size = 18),
                legend.title = element_text(size = 22),
                plot.title = element_text(hjust = 0.5),
                #legend.position = "bottom",
                #legend.spacing.x = unit(0.3, 'cm'),
                title = element_text(size=18),
                plot.margin=grid::unit(c(3,2,3,2), "mm"))+
          #geom_point(shape = 1,size = 3,colour = "black")+
          ggtitle(label = ifelse(length(k) == 1, paste("k =",k), gene.name))
      } else {
        Plots <- list()
        for(k.ind in 1:length(k)){
          Plots[[k.ind]] <- local({
            x.bar <- colMeans(x$x[x$Cs == k[k.ind], which(x$Ds %in% r)])
            ggplot(Coord[-which(x$Ds %in% r),], aes(x, y))+
              theme_bw()+
              geom_point(data = Coord[-which(x$Ds %in% r),], mapping = aes(x, y), size = 3, fill = "white", colour = "gray74", shape = 21)+
              geom_point(data = Coord[which(x$Ds %in% r),], mapping = aes(x, y, fill = x.bar), colour = "white", size = 4, shape = 21)+
              scale_fill_distiller(type = "seq", palette = "Spectral", direction = -1)+
              labs(col = "")+
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.text.x = element_blank(),#element_text(size=18),
                    axis.text.y = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.title=element_blank(),#element_text(size=18),
                    legend.text = element_text(size = 18),
                    legend.title = element_text(size = 22),
                    plot.title = element_text(hjust = 0.5),
                    title = element_text(size=18),
                    plot.margin=grid::unit(c(3,2,3,2), "mm"))+
              #geom_point(shape = 1,size = 3,colour = "black")+
              ggtitle(label = paste("k =",k[k.ind]))})
        }
      }
    }
  }

  Plots
}

