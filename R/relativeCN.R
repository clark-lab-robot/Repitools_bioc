setGeneric("relativeCN", function(input.windows, input.counts, ...)
                                 {standardGeneric("relativeCN")})

setMethod("relativeCN", c("GRanges", "matrix"),
    function(input.windows, input.counts, gc.params = NULL, regions = input.windows[NULL],
              ..., verbose = TRUE)
{
    if(length(input.windows) != nrow(input.counts))
        stop("Rows of counts differ to rows of input regions.\n")

    if(!class(gc.params) %in% c("GCAdjustParams", "NULL"))
        stop("'gc.params' is not an object of class 'GCAdjustParams' or NULL.")
    
    require(GenomicRanges)
    require(DNAcopy)

    if(!is.null(gc.params)) # Do mappability / GC bias adjustment on counts.
    {
        adj.CN <- GCadjustCopy(input.windows, input.counts, gc.params,
                                     verbose = verbose)
        input.scores <- adj.CN@cn
    } else { 
        input.scores <- input.counts
    }

    usable.windows <- !is.na(rowSums(input.scores))
    usable.scores <- input.scores[usable.windows, ]
    usable.locs <- input.windows[usable.windows]

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
    if(verbose == TRUE) message("Segmenting and smoothing copy number ratio.")
    cn <- segment(smooth.CNA(cn), weights = wts, ..., verbose = 0)
    if(verbose == TRUE) message("Done segmenting and smoothing.")
    # Extend CNV region to the end of the interval, since all positions are starts.
    cn$out[, "loc.end"] <- cn$out[, "loc.end"] + width(input.windows)[1]
    CNV.windows <- GRanges(cn$out[, "chrom"],
                           IRanges(cn$out[, "loc.start"], cn$out[, "loc.end"]))

    if(verbose) message("Mapping copy number segmentation to input windows.")
    map <- findOverlaps(input.windows, CNV.windows, select = "first")
    relative.cn <- 2^cn$out[map, "seg.mean"]

    if(length(regions) > 0)
    {	
        if(verbose) message("Mapping copy number windows to regions.")
        map <- findOverlaps(regions, input.windows, select = "first")
        if(is.null(gc.params))
        {
            relative.cn <- relative.cn[map]
        } else {
            adj.CN@windows <- regions
            adj.CN@raw.counts <- adj.CN@raw.counts[map, , drop = FALSE]
            adj.CN@mappability <- adj.CN@mappability[map]
            adj.CN@gc <- adj.CN@gc[map]
            adj.CN@seg.cn <- matrix(relative.cn[map])
        }
    } else
    {
        if(!is.null(gc.params))
            fc.matrix <- matrix(relative.cn)
            colnames(fc.matrix) <- "Fold Change"
            adj.CN@seg.cn <- fc.matrix
    }

    if(is.null(gc.params)) relative.cn else adj.CN
})

setMethod("relativeCN", c("data.frame", "matrix"),
    function(input.windows, input.counts, gc.params = NULL, regions = input.windows[NULL, ],
              ..., verbose = TRUE)
{
    relativeCN(annoDF2GR(input.windows), input.counts, gc.params, annoDF2GR(regions),
               ..., verbose = verbose)
})
