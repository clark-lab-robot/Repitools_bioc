setGeneric("plotClusters", function(x, s.col, non.cl, ...)
          {standardGeneric("plotClusters")})

setMethod("plotClusters", "GRanges", function(x, s.col = NULL, non.cl = NULL, ...)
{
    if(is.null(s.col))
        stop("Score column not given.")
    if(is.null(non.cl))
        stop("Cluster exclusion ID given.")

    require(GenomicRanges)

    elementMetadata(x) <- DataFrame(elementMetadata(x),
                          TSS = as.numeric(ifelse(strand(x) == '+',
                                           start(x),
                                           end(x))))
    
    x.cl <- split(x, elementMetadata(x)[, "cluster"])
    x.cl <- x.cl[names(x.cl) != non.cl]

    invisible(lapply(x.cl, function(y)
    {
        g.anno <- elementMetadata(y)
        par(oma = c(3, 2, 6, 2), mar = c(5, 4, 5, 1))
        plot(g.anno$TSS, g.anno[, s.col], type = 'h', xaxt = 'n', xlab = "",
             ylab = colnames(g.anno)[s.col], ...)
        title(seqnames(y)[1], outer = TRUE)
        mtext("gene", side = 1, line = 5)
        axis(1, g.anno$TSS, g.anno$name, las = 2)
        axis(3, g.anno$TSS, TRUE, las = 2)
    }))
})

setMethod("plotClusters", "data.frame", function(x, s.col = NULL, non.cl = NULL, ...)
{
    if(is.null(s.col))
        stop("Score column not given.")
    if(is.null(non.cl))
        stop("Cluster exclusion ID given.")

    require(GenomicRanges)

    s.name <- colnames(x)[s.col]
    summaryGR <- annoDF2GR(x)
    s.col <- which(colnames(elementMetadata(summaryGR)) == s.name)
    plotClusters(summaryGR, s.col, non.cl, ...)
})
