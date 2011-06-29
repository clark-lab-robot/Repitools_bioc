setGeneric("absoluteCN", function(input.windows, input.counts, gc.params, ...)
                                 {standardGeneric("absoluteCN")})

setMethod("absoluteCN", c("GRanges", "matrix", "GCAdjustParams"),
    function(input.windows, input.counts, gc.params, segment.sqrt = TRUE,
            regions = input.windows[NULL], ..., verbose = TRUE)
{
    require(GenomicRanges)
    require(DNAcopy)

    adj.CN <- GCadjustCopy(input.windows, input.counts, gc.params, verbose = verbose)
    if(segment.sqrt) adj.CN@cn <- sqrt(adj.CN@cn)

    # Do segmentation.
    if(verbose) message("Segmenting absolute copy number estimates and matching to input windows.")
    adj.CN@cn <- apply(adj.CN@cn, 2, function(x)
                 {
                     cn <- CNA(chrom = as.character(seqnames(input.windows)),
                               maploc = as.numeric(start(input.windows)),
                               genomdat = x)
                     cn <- segment(smooth.CNA(cn), ..., verbose = 0)
                     cn$out[, "loc.end"] <- cn$out[, "loc.end"] + width(input.windows)[1]
                     CNV.windows <- GRanges(cn$out[, "chrom"],
                                            IRanges(cn$out[, "loc.start"], cn$out[, "loc.end"]))
                    map <- findOverlaps(input.windows, CNV.windows, select = "first")			   
                    cn$out[map, "seg.mean"]   
                 })
    if(verbose) message("Done segmenting and matching.")
    
    if(segment.sqrt) adj.CN@cn <- adj.CN@cn^2

    if(length(regions) > 0)
    {
        if(verbose) 
            message("Mapping copy number windows to regions.")
        map <- findOverlaps(regions, input.windows, select = "first")
        
        adj.CN@old.counts <- adj.CN@old.counts[map, ]
        adj.CN@cn <- adj.CN@cn[map, ]
        adj.CN@windows <- regions
        adj.CN@gc <- adj.CN@gc[map]
        adj.CN@mappability <- adj.CN@mappability[map]
    }

    adj.CN
})

setMethod("absoluteCN", c("data.frame", "matrix", "GCAdjustParams"),
    function(input.windows, input.counts, gc.params, regions = input.windows[NULL, ], ...)
{
    absoluteCN(annoDF2GR(input.windows), input.counts, gc.params, regions = annoDF2GR(regions), ...)
})
