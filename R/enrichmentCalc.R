setGeneric("enrichmentCalc", function(x, ...) standardGeneric("enrichmentCalc"))

setMethod("enrichmentCalc", "GRangesList",
    function(x, verbose = TRUE, ...)
{
    samp.names <- if(is.null(names(x))) 1:length(x) else names(x)
    ans <- lapply(1:length(x), function(i)
        {
            if(verbose)
                message("Calculating enrichment in ", samp.names[i], '.')
            enrichmentCalc(x[[i]], verbose, ...)
        })
    ans
})

setMethod("enrichmentCalc", "GRanges",
    function(x, seq.len = NULL, verbose = TRUE)
{
    require(GenomicRanges)

    if(any(is.na(seqlengths(x))))
        stop("Some chromosome lengths missing in Seqinfo of reads.")
    if(!is.null(seq.len))
    {
        if(verbose) message("Resizing sample reads to fragment length.")
        x <- suppressWarnings(resize(x, seq.len))
    }
    if(verbose) message("Getting coverage.")
    cov.table <- colSums(table(coverage(x)))
    if(verbose) message("Tabulating coverage.")
    cov.table <- data.frame(as.numeric(names(cov.table)), cov.table)
    colnames(cov.table) <- c("coverage", "bases")
    cov.table
})
