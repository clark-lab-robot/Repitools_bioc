setGeneric("enrichmentPlot", function(x, ...) standardGeneric("enrichmentPlot"))

setMethod("enrichmentPlot", "GRangesList",
    function(x, seq.len, cols = rainbow(length(x)), xlim = c(0, 20), main = "Enrichment Plot",
             total.lib.size = TRUE, verbose = TRUE, ...)
{
    if (length(cols) != length(x)) stop("'x' and 'cols' must have the same number of elements.")
    if (verbose) message("Calculating enrichment.")
    x.enrich <- enrichmentCalc(x, seq.len, verbose)
    if (total.lib.size) {
        if (verbose) message("Normalising to reads per lane.")
        x.counts <- elementLengths(x)
        for (i in 1:length(x))
            x.enrich[[i]]$coverage <- x.enrich[[i]]$coverage/(x.counts[[i]]/1000000)
    }
    x.text <- "Enrichment Level of reads"
    if(total.lib.size) x.text <- paste("Normalised", x.text)
    plot(x = x.enrich[[1]]$coverage, y = x.enrich[[1]]$bases, type = "l", col = cols[1],
         xlim = xlim, main = main, ylab = "Frequency", log = "y", xlab = x.text, ...)
    if (length(x) > 1) for (i in 2:length(x)) {
	    lines(x=x.enrich[[i]]$coverage, y = x.enrich[[i]]$bases, col = cols[i], ...)
    }
    legend("topright", lty = 1, col = cols, legend = names(x), ...)
    invisible(x.enrich)
})
