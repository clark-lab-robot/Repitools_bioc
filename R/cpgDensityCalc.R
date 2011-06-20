setGeneric("cpgDensityCalc", function(x, organism, ...) standardGeneric("cpgDensityCalc"))

setMethod("cpgDensityCalc", c("GRangesList", "BSgenome"),
    function(x, organism, verbose = TRUE, ...)
{
    samp.names <- if (is.null(names(x))) 1:length(x) else names(x)
    ans <- lapply(1:length(x), function(i)
        {
            if (verbose) message("Getting CpG density for ", samp.names[i])
            cpgDensityCalc(x[[i]], organism, verbose = verbose, ...)
        })
    ans
})

setMethod("cpgDensityCalc", c("GRanges", "BSgenome"),
    function(x, organism, seq.len = NULL, window = NULL,
             w.function = c("none", "linear", "exp", "log"), verbose = TRUE)
{
    require(GenomicRanges)

    w.function <- match.arg(w.function)
    if(!is.null(seq.len)) x <- suppressWarnings(resize(x, seq.len))
    if(!is.null(window)) x <- suppressWarnings(resize(x, window, fix = "center"))
    if(w.function == "none") {
        cpgDensity <- sequenceCalc(x, organism, DNAString("CG"))
    } else {
        CGfinds <- sequenceCalc(x, organism, DNAString("CG"), positions = TRUE)
        CGfinds <- lapply(CGfinds, function(u) if(!is.null(u)) abs(u-window/2) else u)
        if(w.function == "linear") {
            cpgDensity <- sapply(CGfinds, function(d) sum(1-(d/(window/2))))
        } else if(w.function == "log") {
            cpgDensity <- sapply(CGfinds, function(d) sum(log2(2-(d/(window/2)))))
        } else {
            cpgDensity <- sapply(CGfinds, function(d) sum(exp(-5*d/(window/2))))	
        }
        rm(CGfinds)
    }    
    if(verbose) message("CpG density calculated for a sample.")
    return(cpgDensity)
})

setMethod("cpgDensityCalc", c("data.frame", "BSgenome"),
    function(x, organism, ...)
{
    require(GenomicRanges)

    if (is.null(x$position)) x$position <- ifelse(x$strand == '+', x$start, x$end)
    x <- GRanges(x$chr, IRanges(x$position, width=1),
                 seqlengths=seqlengths(organism)[unique(x$chr)])
    cpgDensityCalc(x, organism=organism, ...)
})
