setGeneric("mappabilityCalc", function(x, ...){standardGeneric("mappabilityCalc")})
    
setMethod("mappabilityCalc", "GRanges", function(x, organism, verbose = TRUE)
{
    require(GenomicRanges)

    if(verbose) message("Calculating mappability.")
    strand(x) <- "+"
    chrs <- as.character(seqnames(x))
    regions.by.chr <- split(x, chrs)
    
    mappability.by.chr <- lapply(regions.by.chr, function(y)
    {
        # Handle case of windows overlapping past ends of chromosome.
        chr.name <- as.character(seqnames(y))[1]
        chr.max <- length(organism[[chr.name]])
        which.outside.start <- end(y) < 1
        which.outside.end <- start(y) > chr.max
        which.in.chr <- !which.outside.start & !which.outside.end
        which.past.start <- start(y) < 1 & which.in.chr
        which.past.end <- end(y) > chr.max & which.in.chr
        start.Ns <- -start(y)[which.past.start] + 1
        end.Ns <- end(y)[which.past.end] - chr.max
        start(y)[which.past.start] = 1
        end(y)[which.past.end] = chr.max
        temp <- character()
        temp[which.in.chr] <- getSeq(organism, y[which.in.chr], as.character = TRUE)
        temp[!which.in.chr] <- 'N'
        temp[which.past.start] <- unlist(mapply(function(y, z)
                                  {
                                     paste(paste(rep('N', z), collapse = ''), y, sep = '')
                                  }, temp[which.past.start], start.Ns))
        temp[which.past.end] <- unlist(mapply(function(y, z)
                                {
                                    paste(y, paste(rep('N', z), collapse = ''), sep = '')
                                }, temp[which.past.end], end.Ns))
        
        tempAlphabet <- alphabetFrequency(DNAStringSet(temp), as.prob = TRUE)
        1 - tempAlphabet[, "N"]
    })
    unsplit(mappability.by.chr, chrs)
})
    
setMethod("mappabilityCalc", "data.frame", function(x, organism, window = NULL, ...)
{
    require(GenomicRanges)

    if (is.null(x$position)) x$position <- ifelse(x$strand == '+', x$start, x$end)
    x <- GRanges(x$chr, IRanges(x$position, width = 1))
    if(!is.null(window)) x <- resize(x, window, "center")
    mappabilityCalc(x, organism, ...)
})
