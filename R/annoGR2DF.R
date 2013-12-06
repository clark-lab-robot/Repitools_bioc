setGeneric("annoGR2DF", signature = "anno", function(anno, ...)
                                {standardGeneric("annoGR2DF")})

setMethod("annoGR2DF", "GRanges", function(anno)
{
    annoDF <- as.data.frame(anno)
    colnames(annoDF)[1] <- "chr"
    if('*' %in% annoDF$strand)
        annoDF <- annoDF[, -match("strand", colnames(annoDF))]

    annoDF
})

