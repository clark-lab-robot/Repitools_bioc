setGeneric("relativeCN", function(input.windows, input.counts, ...)
                                 {standardGeneric("relativeCN")})

setMethod("relativeCN", c("GRanges", "matrix"),
    function(input.windows, input.counts, gc.params = NULL, ..., verbose = TRUE)
{
    if(length(input.windows) != nrow(input.counts))
        stop("Rows of counts differ to rows of input regions.\n")

    if(!class(gc.params) %in% c("GCAdjustParams", "NULL"))
        stop("'gc.params' is not an object of class 'GCAdjustParams' or 'NULL'.")
    
    require(GenomicRanges)
    require(DNAcopy)

    if(!is.null(gc.params)) # Do mappability / GC bias adjustment on counts.
    {
        CN.result <- GCadjustCopy(input.windows, input.counts, gc.params,
                                     verbose = verbose)
        input.scores <- CN.result@adj.CN
    } else { 
        input.scores <- input.counts
    }

    usable.windows <- !is.na(rowSums(input.scores))
    usable.scores <- input.scores[usable.windows, ]
    usable.locs <- input.windows[usable.windows]

    if(is.null(gc.params))
        Mvalues <- log2((input.counts[, 2] / sum(input.counts[, 2])) /
                        (input.counts[, 1] / sum(input.counts[, 1])))
    else
        Mvalues <- log2((usable.scores[, 2]) / (usable.scores[, 1]))

    cn <- CNA(chrom = as.character(seqnames(usable.locs)),
             maploc = as.numeric(start(usable.locs)),
          data.type = "logratio",
           genomdat = Mvalues,
           sampleid = "Fold Change")
    totals <- colSums(usable.scores)
    wts <- ((totals[2] - usable.scores[, 2]) / (usable.scores[, 2] * totals[2])
          + (totals[1] - usable.scores[, 1]) / (usable.scores[, 1] * totals[1]))^-1
    non0 <- which(wts > 0)
    cn <- cn[non0, ]
    wts <- wts[non0]
    if(verbose == TRUE) message("Smoothing and segmenting copy number ratio.")
    cn <- segment(smooth.CNA(cn), weights = wts, ..., verbose = 0)
    # Extend CNV region to the end of the interval, since all positions are starts.
    cn$out[, "loc.end"] <- cn$out[, "loc.end"] + width(input.windows)[1] - 1
    CNV.windows <- GRanges(cn$out[, "chrom"],
                           IRanges(cn$out[, "loc.start"], cn$out[, "loc.end"]),
                           `relative CN` = 2^cn$out[, "seg.mean"])

    if(!is.null(gc.params))
        CN.result@seg.CN <- GRangesList(CNV.windows)

    if(is.null(gc.params)) CNV.windows else CN.result
})

setMethod("relativeCN", c("data.frame", "matrix"),
    function(input.windows, input.counts, gc.params = NULL, ..., verbose = TRUE)
{
    relativeCN(annoDF2GR(input.windows), input.counts, gc.params, ..., verbose = verbose)
})
