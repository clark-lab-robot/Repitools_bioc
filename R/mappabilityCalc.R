setGeneric("mappabilityCalc", function(x, ...){standardGeneric("mappabilityCalc")})
    
setMethod("mappabilityCalc", "GRanges", function(x, organism, window = NULL,
          verbose = TRUE)
{
    require(GenomicRanges)

    if(verbose) message("Calculating mappability.")

    if(!is.null(window))
    {
        x.posns <- IRanges(as.numeric(ifelse(strand(x) == '+', start(x), end(x))),
                           width = 1)
        ranges(x) <- x.posns
        x <- resize(x, window, "center")
    }
    strand(x) <- "+"
    chrs <- as.character(seqnames(x))
    regions.by.chr <- split(x, chrs)
    
    mappability.by.chr <- lapply(regions.by.chr, function(y)
    {
        # Handle case of windows overlapping past ends of chromosome.
        chr.name <- as.character(seqnames(y))[1]
        chr.max <- length(organism[[chr.name]])
        inside.regions <- restrict(y, 1, chr.max, keep.all.ranges = TRUE)
        window.seqs <- getSeq(organism, inside.regions)
        unmap.counts <- alphabetFrequency(DNAStringSet(window.seqs))[, 'N']
        unmap.counts <- unmap.counts + width(y) - width(inside.regions)
        1 - (unmap.counts / width(y))
    })
    unsplit(mappability.by.chr, chrs)
})
    
setMethod("mappabilityCalc", "data.frame", function(x, organism, ...)
{
    require(GenomicRanges)

    if (is.null(x$position)) x$position <- ifelse(x$strand == '+', x$start, x$end)
    x <- GRanges(x$chr, IRanges(x$position, width = 1))
    mappabilityCalc(x, organism, ...)
})
