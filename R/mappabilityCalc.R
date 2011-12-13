setGeneric("mappabilityCalc", function(x, ...){standardGeneric("mappabilityCalc")})
    
setMethod("mappabilityCalc", "GRanges", function(x, organism, window = NULL,
          type = c("block", "TSS", "center"), verbose = TRUE)
{
    require(GenomicRanges)

    type <- match.arg(type)
    if(type == "block" && !is.null(window))
        stop("'window' is meaningless when region type is \"block\".")

    if(type %in% c("TSS", "center") && is.null(window))
        stop("Using a reference point but window size around it was not specified.")

    info <- "Calculating mappability"
    if(type == "block") info <- paste(info, "in supplied blocks.")
    if(type == "TSS") info <- paste(info, window, "bases around TSSs.")
    if(type == "center") info <- paste(info, window, "bases around feature centres.")
    
    if(verbose) message(info)

    if(type %in% c("TSS", "center"))
    {
        if(type == "TSS")
            x.posns <- IRanges(as.numeric(ifelse(strand(x) == '+', start(x), end(x))),
                           width = 1)
        if(type == "center")
            x.posns <- IRanges(as.integer((start(x) + end(x)) / 2), width = 1)
        
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
    
setMethod("mappabilityCalc", "data.frame", function(x, organism, window = NULL,
          type = c("block", "TSS", "center"), ...)
{
    require(GenomicRanges)

    type <- match.arg(type)
    if(type == "block" && !is.null(window))
        stop("'window' is meaningless when region type is \"block\".")
    
    if(type %in% c("TSS", "center") && is.null(window))
        stop("Using a reference point but window size around it was not specified.")

    if(type == "block")
    {
        x <- GRanges(x$chr, IRanges(x$start, x$end))
    } else {
        if (is.null(x$position))
        {
            if(type == "TSS")
                x$position <- ifelse(x$strand == '+', x$start, x$end)
            if(type == "center")
                x$position <- as.integer((x$start + x$end) / 2)
        }
        x <- GRanges(x$chr, IRanges(x$position, width = 1))
        x <- resize(x, window, "center")
    }

    mappabilityCalc(x, organism, window = NULL, type = "block", ...)
})
