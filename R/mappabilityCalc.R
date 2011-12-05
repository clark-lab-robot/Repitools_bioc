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
    chr.maxs <- seqlengths(organism)[names(regions.by.chr)]
    
    mappability.by.chr <- mapply(function(y, z)
    {
        # Handle case of windows overlapping past ends of chromosome.
        inside.regions <- restrict(y, 1, z, keep.all.ranges = TRUE)
        window.seqs <- suppressWarnings(getSeq(organism, inside.regions))
        unmap.counts <- alphabetFrequency(DNAStringSet(window.seqs))[, 'N']
        unmap.counts <- unmap.counts + width(y) - width(inside.regions)
        1 - (unmap.counts / width(y))
    }, as.list(regions.by.chr), as.list(chr.maxs), SIMPLIFY = FALSE)
    unsplit(mappability.by.chr, chrs)
})
    
setMethod("mappabilityCalc", "data.frame", function(x, organism, ...)
{
    require(GenomicRanges)

    if (is.null(x$position)) x$position <- ifelse(x$strand == '+', x$start, x$end)
    x <- GRanges(x$chr, IRanges(x$position, width = 1))
    mappabilityCalc(x, organism, ...)
})
