setGeneric("clusterPlots", function(c.list, ...){standardGeneric("clusterPlots")})

setMethod("clusterPlots", "ClusteredScoresList",
	function(c.list, plot.ord = 1:length(c.list), plot.type =
                 c("heatmap", "line", "by.cluster"), heat.bg.col = "black",
                 summarize = c("mean", "median"), symm.scale = FALSE, cols = NULL,
                 t.name = NULL, verbose = TRUE, ...)
{
    c.list <- c.list[plot.ord]
    plot.type <- match.arg(plot.type)
    summarize <- match.arg(summarize)

    n.marks <- length(c.list)

    if(is.null(cols) == TRUE)
    {
	require(gplots)
	if(plot.type %in% c("by.cluster", "line"))
	{
	    cols <- colorpanel(n.marks, "blue", "green", "red")
	} else {
	    cols <- colorpanel(100, "blue", "white", "red")
	}
    }

    cvgs <- tables(c.list)

    # Get summary expression for each cluster. Find descending order.
    expr <- c.list@expr
    expr.name <- c.list@expr.name
    cl.id <- c.list@cluster.id
    cl.levels <- levels(factor(cl.id))
    n.clusters <- length(cl.levels)

    keep <- !is.na(cl.id)
    cvgs <- lapply(cvgs, function(x) x[keep, ])
    cl.id <- cl.id[keep]

    if(!is.null(expr))
    {
        expr <- expr[keep]
        if(summarize == "mean")
            cl.expr <- tapply(expr, factor(cl.id), mean, na.rm = TRUE)
        else
            cl.expr <- tapply(expr, factor(cl.id), median, na.rm = TRUE)
        cl.ord <- order(cl.expr, decreasing = TRUE)
    } else {
        cl.expr <- numeric()
        cl.ord <- 1:n.clusters
    }

    # Get x-axis pos and x, score labels.
    pos.labels <- colnames(cvgs[[1]])
    pos <- as.integer(gsub('%', '', pos.labels)) # Get raw position if labels have
                                                 # percentage signs.

    if(symm.scale)
        score.labels <- c("-100%", "0%", "100%")
    else
        score.labels <- c("0%", "50%", "100%")
    
    if(verbose) message("Generating plot.")
    if(plot.type == "by.cluster")
    {
	# Group each cluster from all epigenetic marks.

	profiles <- list() 
	for(i in 1:n.clusters)
	    profiles[[i]] <- sapply(cvgs, function(x)
                                    if(summarize == "mean")
                                        colMeans(x[cl.id == cl.levels[cl.ord[i]], , drop = FALSE], na.rm = TRUE)
                                    else
                                        apply(x[cl.id == cl.levels[cl.ord[i]], , drop = FALSE], 2, median, na.rm = TRUE)
                                   )

	# Plot the lineplots by cluster.
	invisible(mapply(function(x, y)
	{
            y.max = max(unlist(x)) * 1.1
            if(symm.scale) y.min = -y.max else y.min = 0
            layout(rbind(2:1), widths = c(3, 1))
           
	    par(mai = c(1, 0, 0.8, 0))
	    plot.new()
            legend("topleft", legend = names(c.list), title = "Mark", col = cols, ...)
	    par(mai = c(1, 1, 0.8, 0.1))
            if(!is.na(cl.expr[y]))
            {
                cl.text <- paste('(',
                                ifelse(summarize == "mean", "Mean", "Median"),
                                " Expression: ", round(cl.expr[y], 2), ')', sep = '')
            } else {
                cl.text <- paste("Cluster", cl.levels[y])
            }
	    matplot(pos, x, ylim = c(y.min, y.max), type = 'l', col = cols,
                    xlab = "Relative Position", ylab = "Score", yaxt = 'n',
                    xaxs = 'i', yaxs = 'i', ...)
	    axis(2, at = c(y.min, (y.min + y.max) / 2, y.max), label = score.labels)

            mtext(paste("Within Cluster Coverage", cl.text), outer = TRUE, line = -2)
	}, profiles, as.list(cl.ord)))
    } else if(plot.type == "line") # Plot a table of lineplots.
    {
        old.par <- par(no.readonly = TRUE)
        par(oma = c(5, 6, 5, 2), mar = c(0, 0, 0, 0), xpd = NA)
        if(is.null(c.list@.old.ranges))
            ranges <- lapply(cvgs, range, na.rm = TRUE)
        else
            ranges <- c.list@.old.ranges

	      profiles <- lapply(cvgs, function(x)
                           {
                               lapply(cl.ord, function(y)
                               {
                                   if(summarize == "mean")
                                       colMeans(x[cl.id == y, , drop = FALSE], na.rm = TRUE)
                                   else
                                       apply(x[cl.id == y, , drop = FALSE], 2, function(z) median(z), na.rm = TRUE)
                               })
                           })

        n.boxes <- (n.marks) * n.clusters
        layout(matrix(c(1:n.boxes, rep(n.boxes + 1, n.clusters), if(!is.null(expr)) rep(n.boxes + 2, n.clusters)),
                      nrow = n.clusters),
               widths = c(rep(0.8 / n.marks, n.marks), 0.05, if(!is.null(expr)) 0.15))

        x.lims <- c(pos[1], pos[length(pos)])
        for(col.index in 1:n.marks)
        {
            for(row.index in 1:n.clusters)
            {
                matplot(pos, profiles[[col.index]][[row.index]], type = 'l', ylim = ranges[[col.index]],
                        xlab = "", ylab = "", xaxs = 'i', yaxs = 'i', xaxt = 'n', yaxt = 'n',
                        col = cols[col.index], ...)
                if(row.index == n.clusters && col.index == 1)
                {
                    axis(1, pos, pos.labels)
                    mtext("Relative Position", 1, 3)
                }
                if(row.index == 1 && col.index == 1)
                {
                    labels.ys <- par("usr")[3:4]
                    axis(2, c(labels.ys[1], mean(labels.ys), labels.ys[2]), labels = score.labels, las = 2)
                    mtext("Score", 2, 3)
                }
                if(row.index == 1)
                    title(names(c.list)[col.index], line = 2)
            }    
        }

        plot.new()
        plot.window(xlim = c(0, 2), ylim = c(0, n.clusters), xaxs = 'i', yaxs = 'i')
        title("ID", line = 2)
        text(rep(1, n.clusters), (n.clusters:1) - 0.5, 1:n.clusters)

        if(!is.null(expr))
        {
            box.data <- lapply(cl.levels[rev(cl.ord)], function(x) expr[cl.id == x])
            box.data <- boxplot(box.data, plot = FALSE)

            plot.locs <- c(unlist(box.data$stats), unlist(box.data$out))
            max.expr <- 1.1 * max(plot.locs)
            min.data <- min(plot.locs)
            min.expr <- ifelse(min.data < 0, min.data * 1.1, min.data * 0.9)
            
            plot.new()
            plot.window(xlim = c(min.expr, max.expr), ylim = c(0, n.clusters),
                        xaxs = 'i', yaxs = 'i')
            axis(1)
            if(!is.null(expr.name)) mtext(expr.name, 1, 3)
            bxp(box.data, at = (1:n.clusters) - 0.5, add = TRUE, horizontal = TRUE,
               yaxt = 'n')
        }
        mtext(t.name, font = 2, line = 3, outer = TRUE)
        par(old.par)
    } else { # Plot a heatmap.
        if(is.null(c.list@.old.ranges))
            ranges <- lapply(cvgs, range, na.rm = TRUE)
        else
            ranges <- c.list@.old.ranges

        plot.ranges <- lapply(ranges, function(x)
                                      {
                                         if(symm.scale)
                                             c(-1, 1) * max(abs(x))
                                         else
                                             range(x)
                                      })
	
        par(oma = c(1, 1, 3, 1))
	sort.data <- c.list@sort.data
	sort.name <- c.list@sort.name
	# Get order of all features next.
	if(length(sort.data) == 0)
	    ord <- order(factor(cl.id, levels = cl.levels[rev(cl.ord)]))
	else
	    ord <- order(factor(cl.id, levels = cl.levels[rev(cl.ord)]), sort.data)
	
	# Re-arrange the ChIP and expression data and vector of sort data.
	cvgs <- lapply(cvgs, function(x) x[ord, ])
	if(!is.null(expr)) expr <- expr[ord]
	if(length(sort.data) > 0) sort.data <- sort.data[ord]

	# Plot heatmap.
        extras <- sum(!is.null(expr), !is.null(sort.data))
	switch(as.character(extras),
            `0` = layout(rbind(1:(n.marks + 2)), widths=c(1, rep(3, n.marks)), 0.4),
	    `1` = layout(rbind(1:(n.marks + 3)), widths=c(1, rep(3, n.marks), 0.4, 2)),
	    `2` = layout(rbind(1:(n.marks + 4)), widths=c(1, rep(3, n.marks), 0.4, 2, 1)))
	par(mai = c(1.02, 0.50, 0.82, 0.05))
  	
	n.bins = length(cols)
	image(y=seq(1/n.bins/2, 1-(1/n.bins/2), 1/n.bins), z=rbind(1:n.bins),
              col = cols, axes = FALSE, xlab = "Score", ylab = NA)
	axis(2, at=c(0, 0.5, 1), labels = score.labels, las = 2)

	par(mai=c(1.02,0.05,0.82,0.05))
        cl.sizes <- table(cl.id)[rev(cl.ord)]
	bounds <- cumsum(cl.sizes)
        mapply(function(x, y, z)
	{
	    image(pos, 1:nrow(x), t(x), zlim = z, xlab = "Relative Position",
                  xaxt = 'n', yaxt = 'n', col = cols, main = y)
            usr <- par("usr")
            rect(usr[1], usr[3], usr[2], usr[4], col = heat.bg.col)
	    image(pos, 1:nrow(x), t(x), zlim = z, xaxt = 'n', yaxt = 'n', col = cols, add = TRUE)
            axis(1, pos, labels = pos.labels)		

	    # Add lines delimiting the cluster boundaries.
	    abline(h = bounds[-n.clusters], lwd = 3)
	}, cvgs, names(c.list), plot.ranges)

        cl.midpts <- bounds - cl.sizes / 2
        plot.new()
        plot.window(xlim = c(0, 2), ylim = c(0, bounds[length(bounds)]), xaxs = 'i', yaxs = 'i')
        text(rep(1, n.clusters), cl.midpts, labels = cl.levels[rev(cl.ord)], cex = 2)
        title("ID")

	par(mai = c(1.02, 0.05, 0.82, 0.50))
        if(!is.null(expr))
	    plot(expr, y = 1:length(expr), yaxs = 'i', xlab = expr.name, ylab = NA,
                 yaxt = 'n', ...)
	if(!is.null(sort.data)) plot(sort.data, y = 1:length(sort.data), yaxs = 'i',
                     xlab = sort.name, ylab = NA, yaxt = 'n', ...)	

        if(!is.null(t.name))
	    mtext(t.name, line = 0, font = 2, cex = 1.5, outer = TRUE)
    }
})

setMethod("clusterPlots", "ScoresList", function(c.list, scale = function(x) x,
          cap.q = 0.95, cap.type = c("sep", "all"), n.clusters = NULL,
          plot.ord = 1:length(c.list), expr = NULL, expr.name = NULL,
          sort.data = NULL, sort.name = NULL, plot.type = c("heatmap", "line", "by.cluster"),
          summarize = c("mean", "median"), cols = NULL, t.name = NULL,
          verbose = TRUE, ...)
{
    require(cluster)
    
    c.list <- c.list[plot.ord]
    plot.type <- match.arg(plot.type)
    cap.type <- match.arg(cap.type)
    summarize <- match.arg(summarize)
    cvgs <- tables(c.list)

    if(is.null(n.clusters))
	stop("Number of clusters not given.\n")
    if(!is.null(expr) && nrow(cvgs[[1]]) != length(expr))
	stop("Number of features does not match length of expression vector.\n")
    if(!is.null(sort.data) && is.null(sort.name))
    	stop("'sort.data' provided, 'sort.name' is not.")
    if(!is.null(sort.data) && !is.null(expr) && length(expr) != length(sort.data))
	stop("'sort.data' length not the same as number of features in coverage",
             "matrices or expression data.\n")

    # Get rid of any rows that have 1/2 NAs or more.
    rows.na <- apply(cvgs[[1]], 1, function(x) length(which(is.na(x))))
    drop <- rows.na >= ncol(cvgs[[1]]) / 2
    
    cvgs <- lapply(cvgs, scale)
    if(cap.type == "all") max.cvg <- quantile(abs(do.call(cbind, cvgs)), cap.q,
                                              na.rm = TRUE)

    # Find the maximum score allowable, then make any scores bigger be the
    # maximum score.
    cvgs <- lapply(cvgs, function(x)
    {
        if(cap.type == "sep") max.cvg <- quantile(abs(x), cap.q, na.rm = TRUE)
        x[x > max.cvg] = max.cvg
        x[x < -max.cvg] = -max.cvg
        x / max.cvg
    })

    # Do the some kind of clustering for all marks together.
    all <- do.call(cbind, cvgs)[!drop, ]
    any.na <- any(is.na(all))
    set.seed(100)
    if(any.na)
    {
        if(verbose) message("Doing PAM clustering.")
        cl.id <- rep(NA, length(c.list@anno))
        cl.id[!drop] <- pam(all, n.clusters, cluster.only = TRUE, do.swap = FALSE)
    } else { # No NAs in the data. Use k-means for speed.
        if(verbose) message("Doing k-means clustering.")
        cl.id <- kmeans(all, n.clusters, iter.max = 100)$cluster
    }

    # Make sure order of IDs correspond to increasing expression, if present.
    if(!is.null(expr))
    {
        if(summarize == "mean")
            cl.expr <- tapply(expr, factor(cl.id), mean, na.rm = TRUE)
        else
            cl.expr <- tapply(expr, factor(cl.id), median, na.rm = TRUE)
        cl.ord <- order(cl.expr, decreasing = TRUE)
        cl.id <- sapply(cl.id, function(x) match(x, cl.ord))
    }

    csl <- ClusteredScoresList(c.list, scores = cvgs, cluster.id = cl.id,
                                 expr = expr, expr.name = expr.name,
                                 sort.data = sort.data, sort.name = sort.name)
    clusterPlots(csl, plot.type = plot.type, summarize = summarize, cols = cols,
                 t.name = t.name, verbose = verbose, ...)
    invisible(csl)
})
