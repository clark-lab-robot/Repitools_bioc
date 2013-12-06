setGeneric("relativeCN", function(input.windows, input.counts, ...)
                                 {standardGeneric("relativeCN")})

setMethod("relativeCN", c("GRanges", "matrix"),
    function(input.windows, input.counts, gc.params = NULL, ..., verbose = TRUE)
{
    if(length(input.windows) != nrow(input.counts))
        stop("Rows of counts differ to rows of input regions.\n")

    if(!class(gc.params) %in% c("GCAdjustParams", "NULL"))
        stop("'gc.params' is not an object of class 'GCAdjustParams' or 'NULL'.")
    if(!is.null(gc.params)) # Do mappability / GC bias adjustment on counts.
    {
        CN.result <- GCadjustCopy(input.windows, input.counts, gc.params,
                                     verbose = verbose)
        input.scores <- CN.result@adj.CN
        usable.windows <- !is.na(rowSums(input.scores))
        usable.scores <- input.scores[usable.windows, ]
        usable.locs <- input.windows[usable.windows]
        Mvalues <- list(log2((input.counts[, 2]) / (input.counts[, 1])),
                       log2((usable.scores[, 2]) / (usable.scores[, 1])))
        loc.list <- list(input.windows, usable.locs)
        scores.list <- list(CN.result@unadj.CN, usable.scores)
    } else {
        Mvalues <- list(log2((input.counts[, 2] / sum(input.counts[, 2])) /
                        (input.counts[, 1] / sum(input.counts[, 1]))))
        loc.list <- list(input.windows)
        scores.list <- list(input.counts)
    }

    segments.list <- mapply(function(x, y, z)
    {
        cn <- CNA(chrom = as.character(seqnames(x)),
                 maploc = as.numeric(start(x)),
              data.type = "logratio",
             genomdat = z,
             sampleid = "Fold Change")
        totals <- colSums(na.omit(y))
        wts <- ((totals[2] - y[, 2]) / (y[, 2] * totals[2])
              + (totals[1] - y[, 1]) / (y[, 1] * totals[1]))^-1
        non0 <- which(wts > 0)
        cn <- cn[non0, ]
        wts <- wts[non0]
        if(verbose == TRUE) message("Smoothing and segmenting copy number ratio.")
        cn <- segment(smooth.CNA(cn), weights = wts, verbose = 0)
        # Extend CNV region to the end of the interval, since all positions are starts.
        cn$out[, "loc.end"] <- cn$out[, "loc.end"] + width(x)[1] - 1
        CNV.windows <- GRanges(cn$out[, "chrom"],
                               IRanges(cn$out[, "loc.start"], cn$out[, "loc.end"]),
                               `relative CN` = 2^cn$out[, "seg.mean"])        
    }, loc.list, scores.list, Mvalues)

    fc.label <- paste(colnames(input.counts)[2], '/', colnames(input.counts)[1])
    if(!is.null(gc.params))
    {
        CN.result@unadj.CN.seg <- GRangesList(segments.list[[1]])
        CN.result@adj.CN.seg <- GRangesList(segments.list[[2]])
        CN.result@unadj.CN <- CN.result@unadj.CN[, 2, drop = FALSE] / CN.result@unadj.CN[, 1, drop = FALSE]
        colnames(CN.result@unadj.CN) <- fc.label
        CN.result@adj.CN <- CN.result@adj.CN[, 2, drop = FALSE] / CN.result@adj.CN[, 1, drop = FALSE]
        colnames(CN.result@adj.CN) <- fc.label
        CN.result@type <- "relative"
        CN.result
    } else {
        fc.matrix <- matrix(2^Mvalues[[1]])
        colnames(fc.matrix) <- fc.label
        CopyEstimate(input.windows, fc.matrix, GRangesList(segments.list[[1]]), type = "relative")
    }
})

setMethod("relativeCN", c("data.frame", "matrix"),
    function(input.windows, input.counts, gc.params = NULL, ..., verbose = TRUE)
{
    relativeCN(annoDF2GR(input.windows), input.counts, gc.params, ..., verbose = verbose)
})
