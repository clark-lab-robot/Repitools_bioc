setGeneric("mappabilityCalc", function(x, ...){standardGeneric("mappabilityCalc")})
    
setMethod("mappabilityCalc", "GRanges", function(x, organism, verbose = TRUE)
{
    require(GenomicRanges)

    if(verbose) message("Calculating mappability.")
    strand(x) <- "+"
    chrs <- as.character(seqnames(x))
    regionsByChr <- split(x, chrs)
    mappabilityByChr <- lapply(regionsByChr, function(y)
    {
        temp <- getSeq(organism, y, as.character = FALSE)
        tempAlphabet <- alphabetFrequency(temp, as.prob = TRUE)
        1 - tempAlphabet[, "N"]
    })
    unsplit(mappabilityByChr, chrs)
})
    
setMethod("mappabilityCalc", "data.frame", function(x, organism, window = NULL, ...)
{
    require(GenomicRanges)

    if (is.null(x$position)) x$position <- ifelse(x$strand == '+', x$start, x$end)
    x <- GRanges(x$chr, IRanges(x$position, width = 1))
    if(!is.null(window)) x <- resize(x, window, "center")
    mappabilityCalc(x, organism, ...)
})
