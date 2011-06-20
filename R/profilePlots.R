setGeneric("profilePlots", signature = "x", function(x, ...){standardGeneric("profilePlots")})

setMethod("profilePlots", "ScoresList",
    function(x, summarize = c("mean", "median"), gene.lists, n.samples = 1000,
             confidence = 0.975, legend.plot = "topleft", cols = rainbow(length(gene.lists)),
             verbose = TRUE, ...)
{
    # Various checks.
    for(i in 1:length(gene.lists))
    {
        if(class(gene.lists[[i]]) == "logical")
        {
            if(length(gene.lists[[i]]) != length(x@anno)) 
                stop("boolean gene.lists ", names(gene.lists)[i], " length must equal",
                     " num of rows in annotation")
        } else if (class(gene.lists[[i]]) == "integer") {
	    if(max(gene.lists[[i]]) > length(x@anno)) 
		stop("Element index in gene.lists",  names(gene.lists)[i], " greater",
                     " than num of rows in annotation")
	} else stop("gene.lists elements must a be boolean or integer vector")
    }
    stopifnot(confidence > 0.5, confidence < 1)
    summarize <- match.arg(summarize)

    scores <- tables(x)
    x.vals <- as.numeric(gsub('%', '', colnames(scores[[1]])))

    if (length(legend.plot) != length(scores))
        if (length(legend.plot) != 1)
            stop("legend.plot must be either same length as columns in x or 1")
        else
            legend.plot <- rep(legend.plot, length(scores))

    samp.size <- max(sapply(gene.lists, function(u)
                                        if(class(u) == "logical") sum(u)
                                        else length(u)))

    invisible(mapply(function(s.mat, s.name, l.pos)
    {
        if(verbose) message("Processing sample ", s.name, '.')
        sets.scores <- lapply(gene.lists, function(u) s.mat[u, ])
        inds <- lapply(seq_len(n.samples), function(u) sample(nrow(s.mat), samp.size))
        rand.scores <- lapply(inds, function(u) s.mat[u, ])
        if (summarize == "mean")
        {
            sets.scores <- sapply(sets.scores, function(u) apply(u, 2, mean, na.rm = TRUE))
            rand.scores <- sapply(rand.scores, function(u) apply(u, 2, mean, na.rm = TRUE))
        } else {
	    sets.scores <- sapply(sets.scores, function(u) apply(u, 2, median, na.rm = TRUE))
            rand.scores <- sapply(rand.scores, function(u) apply(u, 2, median, na.rm = TRUE))
        }

	rand.conf <- apply(rand.scores, 1, quantile, p = c(1-confidence, 0.5, confidence))
	
	# Plotting.
        # Set up coordinates for polygon.
	matplot(x.vals, cbind(t(rand.conf), sets.scores), type = 'n', xlab = "Relative Position",
                ylab = "Signal", main = s.name, ...)
        # Plot ploygon first, because it is not transparent.
	polygon(x = c(x.vals, rev(x.vals)), y = c(rand.conf[1, ], rev(rand.conf[3, ])),
               col = "lightblue")
        # Plot the profile lines over the polygon.
	matplot(x.vals, cbind(t(rand.conf), sets.scores), type = 'l',
                lty = c(2, 1, 2, rep(1, length(gene.lists))),
                lwd = c(1, 3, 1, rep(3, length(gene.lists))), add = TRUE,
                col = c("blue", "blue", "blue", cols))
	if(!is.na(l.pos))
            legend(l.pos, legend = names(gene.lists), col = cols, lwd = 3)
    }, scores, names(x), legend.plot))
})
