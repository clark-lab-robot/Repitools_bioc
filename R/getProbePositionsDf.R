setOldClass("AffymetrixCdfFile")

setGeneric("getProbePositionsDf", function(cdf, ...)
           {standardGeneric("getProbePositionsDf")})

setMethod("getProbePositionsDf", "AffymetrixCdfFile",
    function(cdf, chrs = NULL, ..., verbose = TRUE)
{
    require(aroma.affymetrix)
	
    ind <- getCellIndices(cdf, ..., useNames = FALSE, unlist = TRUE, verbose = verbose)
    acp <- AromaCellPositionFile$byChipType(getChipType(cdf))
    ch <- acp[ind,1,drop=TRUE]
    sp <- acp[ind,2,drop=TRUE]
    if(is.null(chrs)) chrs <- ch else chrs <- chrs[ch]

    data.frame(chr = chrs,
               position = sp,
               index = ind,
               stringsAsFactors = FALSE)
})
