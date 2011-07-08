setGeneric("chromosomeCNplots", function(copy, ...)
                                 {standardGeneric("chromosomeCNplots")})

.doChrPlot <- function(coords, CNs, y.min, y.max, title, y.lab, pch, cex, pch.col)
{
    
    plot(x = NA, y = NA, xlim = c(0, max(coords[, 2])), ylim = c(0, y.max),
         xlab = "Position", ylab = y.lab, main = title)

    points(x = as.numeric(t(cbind(coords[, 1], coords[, 2]))), y = rep(CNs, each = 2),
           pch = pch, cex = cex, col = pch.col)
    abline(h = 1:y.max, lty = 2, col = "grey")
}

.drawSegs <- function(segGR, curr.chr, lty, lwd, seg.col)
{
    seg.coords <- as.data.frame(ranges(segGR))[curr.chr, 1:2]
    seg.vals <- elementMetadata(segGR)[curr.chr, 1]
    apply(cbind(seg.coords, seg.vals), 1, function(segment)
               lines(segment[1:2], rep(segment[3], 2),
               col = seg.col, lty = lty, lwd = lwd))
}

setMethod("chromosomeCNplots", c("CopyEstimate"),
    function(copy, y.max = NULL, pch = 19, cex = 0.2,
             pch.col = "black", seg.col = "red", lty = 1, lwd = 2, verbose = TRUE)
{
    if(is.null(y.max))
        stop("Parameter 'y.max' must be provided.")

    # Only relative fold changes can be unadjusted.
    y.min <- -y.max
    y.lab <- "Fold Change"

    chrs <- factor(as.character(seqnames(copy@windows)))
    coords <- as.data.frame(ranges(copy@windows))[, 1:2]
    unadj.CN <- copy@unadj.CN
    IDs <- colnames(unadj.CN)
    sample.segs <- copy@unadj.CN.seg

    options(scipen = 12)
    invisible(mapply(function(x, y, z)
    {
        str(x)
        str(head(y))
str(z)
        lapply(levels(chrs), function(chr)
        {
            winds.curr.chr <- chrs == chr
            title <- paste(x, chr)
            .doChrPlot(coords[winds.curr.chr, ], y[winds.curr.chr], y.min, y.max, title, y.lab, pch, cex, pch.col)
            .drawSegs(z, as.character(seqnames(z)) == chr, lty, lwd, seg.col)
        })
    }, as.list(IDs),
       split(unadj.CN, rep(1:ncol(unadj.CN), each = nrow(unadj.CN))),
       as.list(sample.segs)))
})

setMethod("chromosomeCNplots", c("AdjustedCopyEstimate"),
    function(copy, y.max = NULL, pch = 19, cex = 0.2,
             pch.col = "black", seg.col = "red", lty = 1, lwd = 2, verbose = TRUE)
{
    if(is.null(y.max))
        stop("Parameter 'y.max' must be provided.")

    if(copy@type == "absolute")
    {
        y.min <- 0 
        y.lab <- "Absolute Copy"
    } else {
        y.min <- -y.max
        y.lab <- "Fold Change"
    }


    chrs <- factor(as.character(seqnames(copy@windows)))
    coords <- as.data.frame(ranges(copy@windows))[, 1:2]
    unadj.CN <- copy@unadj.CN
    IDs <- colnames(unadj.CN)
    unadj.CN <- split(unadj.CN, rep(1:ncol(unadj.CN), each = nrow(unadj.CN)))
    adj.CN <- copy@adj.CN
    adj.CN <- split(adj.CN, rep(1:ncol(adj.CN), each = nrow(adj.CN)))
    unadj.segs <- copy@unadj.CN.seg
    adj.segs <- copy@adj.CN.seg

    options(scipen = 12)
    par(mfrow = 2:1)
    invisible(mapply(function(v, w, x, y, z)
    {
        lapply(levels(chrs), function(chr)
        {
            winds.curr.chr <- chrs == chr
            title <- paste(v, chr, "Before GC Adjustment")
            .doChrPlot(coords[winds.curr.chr, ], w[winds.curr.chr], y.min, y.max, title, y.lab, pch, cex, pch.col)
            .drawSegs(y, as.character(seqnames(y)) == chr, lty, lwd, seg.col)
            title <- paste(v, chr, "After GC Adjustment")
            .doChrPlot(coords[winds.curr.chr, ], x[winds.curr.chr], y.min, y.max, title, y.lab, pch, cex, pch.col)
            .drawSegs(z, as.character(seqnames(z)) == chr, lty, lwd, seg.col)
        })
    }, IDs, unadj.CN, adj.CN, as.list(unadj.segs), as.list(adj.segs)))
})


