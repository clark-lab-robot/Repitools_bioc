setGeneric("relativeCN", function(input.windows, ip.windows, ...){standardGeneric("relativeCN")})

setMethod("relativeCN", c("GRanges", "GRanges"),
    function(input.windows, ip.windows, input.counts = NULL, gc.params = NULL, ..., verbose = TRUE)
{
    if(is.null(input.counts))
        stop("Matrix of counts for input types not given.")

    if(length(input.windows) != nrow(input.counts))
        stop("Rows of counts differ to rows of input regions.\n")
    
    require(GenomicRanges)
    require(DNAcopy)

    if(!is.null(gc.params)) # Do mappability / GC bias adjustment on counts.
        input.scores <- absoluteCN(input.windows, ip.windows = input.windows, input.counts,
                                   gc.params, verbose)
    else 
        input.scores <- input.counts

    usable.windows <- !is.na(rowSums(input.scores))
    usable.scores <- input.scores[usable.windows, ]
    usable.locs <- input.windows[usable.windows, ]

    Mvalues <- log2((usable.scores[, 2] / sum(usable.scores[, 2])) /
                    (usable.scores[, 1] / sum(usable.scores[, 1])))
    cn <- CNA(chrom = as.character(seqnames(usable.locs)),
             maploc = as.numeric(start(usable.locs)),
          data.type = "logratio",
           genomdat = Mvalues,
           sampleid = paste(colnames(usable.scores)[2], "/",
                            colnames(usable.scores)[1], "Fold Change"))
    totals <- colSums(usable.scores)
    wts <- ((totals[2] - usable.scores[, 2]) / (usable.scores[, 2] * totals[2])
          + (totals[1] - usable.scores[, 1]) / (usable.scores[, 1] * totals[1]))^-1
    non0 <- which(wts > 0)
    cn <- cn[non0, ]
    wts <- wts[non0]
    if(verbose == TRUE) message("Segmenting and smoothing.")
    cn <- segment(smooth.CNA(cn), weights = wts,p.method = "perm", undo.splits = "sdundo", undo.SD = 1, verbose = 0)
    if(verbose == TRUE) message("Done segmenting and smoothing.")
    # Extend CNV region to the end of the interval, since all positions are starts.
    cn$out[, "loc.end"] <- cn$out[, "loc.end"] + width(input.windows)[1]
    CNV.windows <- GRanges(cn$out[, "chrom"],
                           IRanges(cn$out[, "loc.start"], cn$out[, "loc.end"]))
	
    if(verbose == TRUE) message("Mapping copy number windows to IP windows.")
    map <- findOverlaps(ip.windows, CNV.windows, select = "first")
							   
    relative.cn <- 2^cn$out[map, "seg.mean"]
    names(relative.cn) <- .getNames(ip.windows)

    relative.cn
})

setMethod("relativeCN", c("data.frame", "data.frame"),
    function(input.windows, ip.windows, input.counts = NULL, gc.params = NULL, ..., verbose = TRUE)
{
    relativeCN(annoDF2GR(input.windows), annoDF2GR(ip.windows), input.counts, gc.params,
               ..., verbose = verbose)
})
