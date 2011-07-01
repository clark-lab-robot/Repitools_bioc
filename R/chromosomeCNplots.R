setGeneric("chromosomeCNplots", function(copy, ...)
                                 {standardGeneric("chromosomeCNplots")})

.doChrPlot <- function(coords, CNs, y.max, title, pch, cex, pch.col)
{
    
    plot(x = NA, y = NA, xlim = c(0, max(coords[, 2])), ylim = c(0, y.max),
         xlab = "Position", ylab = "Copy Estimate", main = title)

    points(x = as.numeric(t(cbind(coords[, 1], coords[, 2]))), y = rep(CNs, each = 2),
           pch = pch, cex = cex, col = pch.col)
    abline(h = 1:y.max, lty = 2, col = "grey")
}

setMethod("chromosomeCNplots", c("CopyEstimate"),
    function(copy, y.max = NULL, pch = 19, cex = 0.2,
             pch.col = "black", seg.col = "red", lty = 1, lwd = 2, verbose = TRUE)
{
    if(is.null(y.max))
        stop("Parameter 'y.max' must be provided.")

    chrs <- factor(as.character(seqnames(copy@windows)))
    coords <- as.data.frame(ranges(copy@windows))[, 1:2]
    unadj.CN <- copy@unadj.CN
    IDs <- colnames(unadj.CN)
    unadj.CN <- split(unadj.CN, rep(1:ncol(unadj.CN), each = nrow(unadj.CN)))
    adj.CN <- copy@adj.CN
    adj.CN <- split(adj.CN, rep(1:ncol(adj.CN), each = nrow(adj.CN)))
    sample.segs <- copy@seg.CN

    options(scipen = 12)
    par(mfrow = 2:1)
    invisible(mapply(function(w, x, y, z)
    {
        lapply(levels(chrs), function(chr)
        {
            curr.chr <- chrs == chr
            segs.chr <- as.character(seqnames(z)) == chr
            title <- paste(w, chr, "Before GC Adjustment")
            .doChrPlot(coords[curr.chr, ], x[curr.chr], y.max, title, pch, cex, pch.col)
            title <- paste(w, chr, "After GC Adjustment")
            .doChrPlot(coords[curr.chr, ], y[curr.chr], y.max, title, pch, cex, pch.col)
            seg.coords <- as.data.frame(ranges(z))[segs.chr, 1:2]
            seg.vals <- elementMetadata(z)[segs.chr, 1]
            apply(cbind(seg.coords, seg.vals), 1, function(segment)
                  lines(segment[1:2], rep(segment[3], 2),
                  col = seg.col, lty = lty, lwd = lwd))
        })
    }, IDs, unadj.CN, adj.CN, as.list(sample.segs)))
})
