setGeneric("annoDF2GR", signature = "anno", function(anno, ...)
                                {standardGeneric("annoDF2GR")})

setMethod("annoDF2GR", "data.frame", function(anno)
{
    require(GenomicRanges)

    col.missing <- setdiff(c("chr", "start", "end"), colnames(anno))
    if(length(col.missing) > 0)
	stop("Columns ", paste(col.missing, collapse = ", "),
            " of annotation are not present.")

    extra <- !colnames(anno) %in% c("chr", "start", "end", "strand")
    DF <- DataFrame(anno[, extra, drop = FALSE])    

    GRanges(anno$chr,
	    IRanges(anno$start, anno$end),
            if("strand" %in% colnames(anno)) anno$strand else '*',
            DF)
})
