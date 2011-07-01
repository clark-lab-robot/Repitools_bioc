setGeneric("GCbiasPlots", function(copy, ...)
                                 {standardGeneric("GCbiasPlots")})

.doBiasPlot <- function(GC, CN, title, y.max, pch, cex, pch.col, lty, lwd, line.col)
{
    plot(GC, CN, ylim = c(0, y.max), xlab = "GC Content", ylab = "Copy Number",
         main = title, pch = pch, cex = cex, col = pch.col)
    lines(lowess(GC, CN), lty = lty, lwd = lwd, col = line.col)    
}

setMethod("GCbiasPlots", c("CopyEstimate"),
    function(copy, y.max = NULL, pch = 19, cex = 0.2,
             pch.col = "black", line.col = "red", lty = 1, lwd = 2, verbose = TRUE)
{
    if(is.null(y.max))
        stop("Parameter 'y.max' must be provided.")

    windows <- copy@windows
    unadj.CN <- copy@unadj.CN
    adj.CN <- copy@adj.CN
    GC <- elementMetadata(windows)[, "GC"]
    mappability <- elementMetadata(windows)[, "mappability"]
    if(length(y.max) == 1 && ncol(unadj.CN) > 1)
        y.max <- rep(y.max, ncol(unadj.CN))

    use <- !is.na(GC) & !is.na(mappability) & rowSums(!is.na(unadj.CN))
    windows <- windows[use]
    unadj.CN <- unadj.CN[use, , drop = FALSE]
    adj.CN <- adj.CN[use, , drop = FALSE]
    GC <- GC[use]
    
    par(mfrow = 2:1)
    par(mar = c(4, 4, 1, 1))

    invisible(mapply(function(sample.unadj.CN, sample.adj.CN, sample.ID, y.max)
    {
        if(verbose) message("Generating GC bias plots for ", sample.ID)
        .doBiasPlot(GC, sample.unadj.CN, paste(sample.ID, "Unadjusted Copy Numbers"), y.max,
                   pch, cex, pch.col, lty, lwd, line.col)
        .doBiasPlot(GC, sample.adj.CN, paste(sample.ID, "GC Adjusted Copy Numbers"), y.max,
                   pch, cex, pch.col, lty, lwd, line.col)
        abline(h = 1:y.max, lty = 2, col = "grey")

    }, split(unadj.CN, rep(1:ncol(unadj.CN), each = nrow(unadj.CN))),
       split(adj.CN, rep(1:ncol(adj.CN), each = nrow(adj.CN))),
       as.list(colnames(unadj.CN)),
       as.list(y.max)))
})
