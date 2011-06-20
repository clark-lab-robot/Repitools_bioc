setGeneric("mergeReplicates", signature = "reads", function(reads, types, ...)
                                         {standardGeneric("mergeReplicates")})

setMethod("mergeReplicates", "GRangesList", function(reads, types, verbose = TRUE)
{
    require(GenomicRanges)

    if(is.null(types))
    	stop("Mandatory argument 'types' not provided.\n")
    if(length(types) != length(reads))
    	stop("'types' and 'reads' lengths differ.\n")

    if(verbose) message("Unlisting GRangesList.\n")
    readsGR <- unlist(reads, use.names = FALSE)
    rdTypes <- Rle(types, elementLengths(reads))
    if(verbose) message("Splitting by types.\n")
    reads <- split(readsGR, rdTypes)
    metadata(reads) <- list(names(reads))
    gc()
    if(verbose) message("Pooled GRangesList created.\n")
    reads
})
