setGeneric("absoluteCN", function(input.windows, input.counts, gc.params, ...)
                                 {standardGeneric("absoluteCN")})

setMethod("absoluteCN", c("GRanges", "matrix", "GCAdjustParams"),
    function(input.windows, input.counts, gc.params, segment.sqrt = TRUE, ...,
             verbose = TRUE)
{
    require(GenomicRanges)
    require(DNAcopy)

    CN.result <- GCadjustCopy(input.windows, input.counts, gc.params, verbose = verbose)
    if(segment.sqrt) CN <- sqrt(CN.result@adj.CN) else CN <- CN.result@adj.CN

    # Do segmentation.
    if(verbose) message("Smoothing and segmenting absolute copy number estimates.")
    CN.result@seg.CN <- GRangesList(apply(CN, 2, function(x)
                        {
                            cna <- CNA(chrom = as.character(seqnames(input.windows)),
                                      maploc = as.numeric(start(input.windows)),
                                      genomdat = x)
                            cna <- segment(smooth.CNA(cna), ..., verbose = 0)
                            if(segment.sqrt) cna$out[, "seg.mean"] <- cna$out[, "seg.mean"]^2
                            cna$out[, "loc.end"] <- cna$out[, "loc.end"] + width(input.windows)[1]
                            GRanges(cna$out[, "chrom"],
                                    IRanges(cna$out[, "loc.start"], cna$out[, "loc.end"]),
                                    CN = cna$out[, "seg.mean"])
                        }))
    CN.result
})

setMethod("absoluteCN", c("data.frame", "matrix", "GCAdjustParams"),
    function(input.windows, input.counts, gc.params, ...)
{
    absoluteCN(annoDF2GR(input.windows), input.counts, gc.params, ...)
})
