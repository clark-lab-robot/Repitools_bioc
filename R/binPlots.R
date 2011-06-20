setGeneric("binPlots", signature = "x", function(x, ...){standardGeneric("binPlots")})

setMethod("binPlots", "ScoresList",
    function(x, summarize = c("mean", "median"), ordering, ord.label,
             plot.type = c("line", "heatmap", "terrain"), n.bins = 10, cols = NULL,
             lwd = 3, lty = 1, same.scale = TRUE, symm.scale = FALSE, verbose = TRUE)
{
    summarize <- match.arg(summarize)
    scores <- tables(x)
    n.points <- ncol(scores[[1]])
    x.vals <- as.numeric(gsub('%', '', colnames(scores[[1]])))
    x.gap <- x.vals[2] - x.vals[1]

    def.par <- par(no.readonly = TRUE) # save default, for resetting.
    plot.type <- match.arg(plot.type)
    if(!ncol(ordering) == length(scores))
    {
        if(!ncol(ordering) == 1)
            stop("Ordering must have either 1 column or the same number as the score tables.")
        ordering <- ordering[, rep(1, length(scores))]
        colnames(ordering) <- gsub("\\.[0-9]+", '', colnames(ordering))
    }

    if(is.null(cols))
    {
        require(gplots)
	if(plot.type == "line") {
	  cols <- colorpanel(n.bins, "blue", "green", "red")
	} else {
	  cols <- colorpanel(64, "blue", "white", "red")
	}
    }
  
    o.types <- character()
    for(i in 1:ncol(ordering))
    {
        if(class(ordering[, i]) == "numeric")
            o.types[i] <- "Order:"
        else if(class(ordering[, i]) == "factor")
            o.types[i] <- "Factor:"
    }
    
    # Split genes into intervals.
    if(verbose) message("Calculating intervals.")
    breaks <- apply(ordering, 2, function(u)
                   {
                       if(class(u) == "numeric") {
                           br <- quantile(u, p = (0:n.bins)/n.bins)
                           list(breakpoints = br, intervals = cut(u, breaks = br))
	               } else if(class(u) == "factor") {
	                   n.bins <- length(levels(u))
	                   list(breakpoints = u, intervals = u)
	               }
                   })

    # Get matrices for each ordering bin, in each dataset.
    if(verbose) message("Splitting scores into bins.")
    items.bins <- mapply(function(y, z)
                         {
                            sub.matrices <- lapply(levels(z$intervals), function(l)
                                                  {
                                                       y[z$intervals == l, ]
                                                  })
                            names(sub.matrices) <- levels(z$intervals)
                            sub.matrices
                         }, scores, breaks, SIMPLIFY = FALSE)

    # Make summaries for plotting.
    items.bins <- lapply(items.bins, function(y)
                        {
                            lapply(y, function(z)
                                  {
                                      if(summarize == "mean")
                                          apply(z, 2, mean, na.rm = TRUE)
                                      else
                                          apply(z, 2, median, na.rm = TRUE)
                                  })
                         })

    # Calculcate ranges for each plot.
    if(same.scale)
        ranges <- rep(list(range(items.bins, na.rm = TRUE)), length(scores))
    else
        ranges <- lapply(items.bins, range, na.rm = TRUE)
    if(symm.scale)
        ranges <- lapply(ranges, function(rng) c(-max(abs(ranges)), max(abs(ranges))))

    if(verbose) message("Drawing plots.")
    invisible(mapply(function(item.bins, item.name, range, o.type, o.name, item.breaks)
    {
        if(verbose) message("Plotting ", item.name, '.')
        cut.levels <- levels(item.breaks[["intervals"]])
        
        t.name <- paste("Signal:", item.name, o.type, o.name, sep = ' ')
        if(plot.type == "line")
        {
            layout(rbind(c(1, 2)), widths = c(3, 2))
            par(oma = c(0, 0, 2, 0))
            par(mai = c(1.02, 0.90, 0.82, 0))
            matplot(x.vals, do.call(cbind, item.bins), type = 'l', col = cols,
                    lty = lty, lwd = lwd, xlab = "Relative Position", ylab = "Signal",
                    ylim = range)
            par(mai = c(1.02, 0.05, 0.82, 0))
            plot.new()
            legend(x = "top", title = ord.label, col = cols, lty = 1, legend = cut.levels)
            mtext(t.name, line = 0.5, outer = TRUE)
        } else if(plot.type == "heatmap") {
               layout(rbind(c(1, 2, 3)), widths=c(1, 3, 1))
               par(mai = c(1.02, 0.50, 0.82, 0.05))
               par(oma = c(0, 0, 0, 0))
               image(y = seq(1/n.bins/2, 1-(1/n.bins/2), 1/n.bins), z = rbind(1:n.bins),
                     col = cols, axes = FALSE, xlab = "Signal Intensity", ylab = NA)
               i.labels <- format(seq(range[1], range[2], length.out = n.bins+1), digits = 1) 
               axis(2, at=(0:n.bins)/n.bins, labels = i.labels)
               par(mai=c(1.02,0.05,0.82,0.05))
               image(x.vals, 1:n.bins, do.call(cbind, item.bins), xlab = "Relative Position",
                     yaxt = 'n', ylab = "Bin", col = cols, zlim = range)
                par(mai = c(1.02, 0.05, 0.82, 0.50))
                bounds <- item.breaks[["breakpoints"]]
                plot(x = bounds, y = 0:n.bins, type = 'l', yaxt = 'n', lwd = 3, xlab = ord.label, yaxs = 'i')
                par(oma = c(0, 0, 2, 0))
		mtext(t.name, line = 0, outer = TRUE)
        } else if(plot.type == "terrain") {
              dm <- do.call(cbind, item.bins)
              layout(1)
              par(oma = c(0, 0, 2, 0))
              dm.avg <- (dm[-1, -1] + dm[-1, -(ncol(dm) - 1)] +
                         dm[-(nrow(dm) -1), -1] + dm[-(nrow(dm) - 1), -(ncol(dm) - 1)]) / 4

              dm.bins <- cut(dm.avg, breaks = seq(range[1], range[2], length.out = length(cols)),
                             include.lowest = TRUE)
              use.cols = cols[dm.bins]
              persp(x.vals, 1:n.bins, dm, xlab = "Relative Position", ylab = ord.label,
                  col = use.cols, zlim = range, theta = -25, phi = 20, d = 1.5, border = NA,
                  ticktype = "detailed", zlab = "Signal")
	      mtext(t.name, line = 0, outer = TRUE)
	}
        par(def.par)#- reset to default
    }, items.bins, names(x), ranges, o.types, colnames(ordering), breaks))
})
