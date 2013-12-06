setGeneric("absoluteCN", function(input.windows, input.counts, gc.params, ...)
                                 {standardGeneric("absoluteCN")})

setMethod("absoluteCN", c("GRanges", "matrix", "GCAdjustParams"),
    function(input.windows, input.counts, gc.params, segment.sqrt = TRUE, ...,
             verbose = TRUE)
{
    CN.result <- GCadjustCopy(input.windows, input.counts, gc.params, verbose = verbose)
    if(segment.sqrt) CN.raw <- sqrt(CN.result@unadj.CN) else CN.raw <- CN.result@unadj.CN
    if(segment.sqrt) CN.adj <- sqrt(CN.result@adj.CN) else CN.adj <- CN.result@adj.CN

    # Do segmentation on both unadjusted and adjusted estimates.
    if(verbose) message("Smoothing and segmenting absolute copy number estimates.")
    segments.list <- lapply(list(CN.raw, CN.adj), function(x)
    {
        GRangesList(apply(x, 2, function(y)
        {
            cna <- CNA(chrom = as.character(seqnames(input.windows)),
                       maploc = as.numeric(start(input.windows)),
                       genomdat = y)
            cna <- segment(smooth.CNA(cna), ..., verbose = 0)
            if(segment.sqrt) cna$out[, "seg.mean"] <- cna$out[, "seg.mean"]^2
            cna$out[, "loc.end"] <- cna$out[, "loc.end"] + width(input.windows)[1]
            GRanges(cna$out[, "chrom"], IRanges(cna$out[, "loc.start"], cna$out[, "loc.end"]),
                    CN = cna$out[, "seg.mean"])
        }))
    })
    CN.result@unadj.CN.seg = segments.list[[1]]
    CN.result@adj.CN.seg = segments.list[[2]]
    CN.result@type = "absolute"
    CN.result
})

setMethod("absoluteCN", c("data.frame", "matrix", "GCAdjustParams"),
    function(input.windows, input.counts, gc.params, ...)
{
    absoluteCN(annoDF2GR(input.windows), input.counts, gc.params, ...)
})
