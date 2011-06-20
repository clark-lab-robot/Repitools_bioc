# Common short functions used by multiple function in Repitools.

setGeneric(".validate", signature = c("anno"), function(anno, up, down)
                                        {standardGeneric(".validate")})

setMethod(".validate", "GRanges", function(anno, up, down)
{
    if(is.null(up))
        stop("Mandatory argument 'up' not given.")
    if(is.null(down))
        stop("Mandatory argument 'down' not given.")

    str <- strand(anno)

    if(-up > down)
	stop("'up' and 'down' window boundaries cross over each other.\n")
	
    if(any(str == '*') && any(str %in% c('+', '-')))
	stop("Annotation contains mixed feature types.")

    # For unstranded features.
    if(any(str %in% '*') && up != down)
	stop("Different upstream and downstream window distances don't make
              sense for unstranded features.\n")
    NULL
})

.getNames <- function(annoGR)
{
    if("name" %in% names(elementMetadata(annoGR)))
        return(elementMetadata(annoGR)[, "name"])
    else
        return(1:length(annoGR))
}
