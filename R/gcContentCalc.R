setGeneric("gcContentCalc", function(x, organism, ...){standardGeneric("gcContentCalc")})

setMethod("gcContentCalc", c("GRanges", "BSgenome"),
    function(x, organism, verbose = TRUE)
{
    require(GenomicRanges)
    require(BSgenome)
    
    if(verbose) message("Calculating GC content.")
    strand(x) <- "+"
    temp <- getSeq(organism, x, as.character=FALSE)
    tempAlphabet <- alphabetFrequency(temp, as.prob=TRUE)
    (tempAlphabet[,"C"]+tempAlphabet[,"G"])
})

setMethod("gcContentCalc", c("data.frame", "BSgenome"),
    function(x, organism, window = NULL, ...)
{
    if(is.null(window))
        stop("Window size not given.")
    require(GenomicRanges)

    if (is.null(x$position)) x$position <- ifelse(x$strand == '+', x$start, x$end)
    x <- GRanges(x$chr, IRanges(x$position, width=1), seqlengths=seqlengths(organism)[unique(x$chr)])
    x <- resize(x, window, fix="center")
    gcContentCalc(x, organism, ...)
})
