# A collection of variables that describe where the sampling will happen.
# An S4 class, so that they are created once, then dispatched on when
# required for multiple samples.

setClass(".CoverageSamples",
         representation(
                        pos.labels = "ANY", # character or numeric.
                        cvg.samps = "GRanges"
			))

setClass("ScoresList", representation(
                                        names = "character",
					scores = "list", # list of matrices.
					anno = "GRanges",
					up = "numeric",
					down = "numeric",
					dist = "ANY", # character or NULL.
					freq = "numeric",
					s.width = "ANY", # character or NULL.
                                        .samp.info = ".CoverageSamples"
					))

setMethod("names", "ScoresList", function(x) x@names)
setGeneric("tables", function(x) {standardGeneric("tables")})
setMethod("tables", "ScoresList", function(x) x@scores)
setMethod("length", "ScoresList", function(x) length(x@scores))

setMethod("show", "ScoresList",
    function(object)
    {
        if(!is.null(object@dist))
	    dist.label <- ifelse(object@dist == "percent", '%', "bases")
        else
            dist.label <- "bases"
	cat("An object of class 'ScoresList'.\n")
	cat("Tables: ", paste(object@names, collapse = ", "), ".\n", sep = '')
	cat("Features:\n")
	print(object@anno)
	cat("Region:",  paste(object@up, dist.label, "up to", object@down,
                         dist.label, "down.\n"))
        if(!is.null(object@s.width))
        {
	    cat("Smoothing:", paste(object@s.width, collapse = ", "), "bases.\n")
	    cat("Sampling : ", object@freq, ' ', dist.label, ".\n\n",  sep = '')
        } else {
            cat("Window Width : ", object@freq, ' ', dist.label, ".\n\n",  sep = '')
        }        
    })

setMethod("[", "ScoresList",
    function(x, i)
    {
	new("ScoresList", names = x@names[i], anno = x@anno, scores = x@scores[i],
	                  up = x@up, down = x@down, dist = x@dist,
			  freq = x@freq, s.width = x@s.width[i], .samp.info = x@.samp.info)
    })

setReplaceMethod("names", "ScoresList",
    function(x, value)
    {
	x@names <- value
	x
    }
)

setGeneric("subsetRows", function(x, i) {standardGeneric("subsetRows")})
setMethod("subsetRows", "ScoresList",
    function(x, i = NULL)
{
    if(is.null(i))
        stop("No row indices given to subset by.")

    new("ScoresList", names = x@names, anno = x@anno[i],
                      scores = lapply(x@scores, function(y) y[i, ]),
                      up = x@up, down = x@down, dist = x@dist,
                      freq = x@freq, s.width = x@s.width, .samp.info = x@.samp.info)
})

setClass("ClusteredScoresList", representation(
                                    cluster.id = "numeric",
				    expr = "ANY",
                                    expr.name = "ANY",
				    sort.data = "ANY",
				    sort.name = "ANY",
                                    .old.ranges = "ANY"),
                       contains = "ScoresList")

setMethod("show", "ClusteredScoresList",
    function(object)
    {
	dist.label <- ifelse(object@dist == "percent", '%', "bases")
	cat("An object of class 'ClusteredScoresList'.\n")
	cat("Tables: ", paste(object@names, collapse = ", "), ".\n", sep = '')
	cat("Region: ",  paste(object@up, dist.label, "up to", object@down,
	    dist.label, "down.\n"))
	cat("Features:\n")
	print(object@anno)
        if(!is.null(object@s.width))
	    cat("Smoothing:", paste(object@s.width, collapse = ", "), "bases.\n")
	cat("Sampling: ", object@freq, ' ', dist.label, ".\n",  sep = '')
        if(!is.null(object@expr))
            cat("Feature Expressions:", object@expr.name,
                paste('\n', paste(head(round(object@expr, 2)), collapse = ", "),
                ", ...\n", sep = ''))
	cat("Feature Clusters:", paste(paste(head(object@cluster.id),
	    collapse = ", "), ", ...\n", sep = ''))
	if(!is.null(object@sort.data))
	    cat("Within Cluster Sorting: By ", object@sort.name, ". ",
	    paste(paste(head(object@sort.data), collapse = ", "), ", ...\n", sep = ''),
	    sep = '')		
    })

# Constructor
setGeneric("ClusteredScoresList", function(x, ...)
           {standardGeneric("ClusteredScoresList")})
setMethod("ClusteredScoresList", "ScoresList",
    function(x, scores = tables(x), expr = NULL, expr.name = NULL, cluster.id, sort.data = NULL,
             sort.name = NULL)
{
	new("ClusteredScoresList", names = x@names, scores = scores, anno = x@anno,
	    up = x@up, down = x@down, dist = x@dist,
	    freq = x@freq, s.width = x@s.width, cluster.id = cluster.id,
	    expr = expr, expr.name = expr.name, sort.data = sort.data,
            sort.name = sort.name, .samp.info = x@.samp.info)
})

setMethod("[", "ClusteredScoresList",
    function(x, i)
{
    new("ClusteredScoresList", names = x@names[i], scores = x@scores[i],
	anno = x@anno, up = x@up, down = x@down, dist = x@dist,
	freq = x@freq, s.width = x@s.width[i], cluster.id = x@cluster.id,
	expr = x@expr, expr.name = x@expr.name, sort.data = x@sort.data,
        sort.name = x@sort.name, .samp.info = x@.samp.info)
    })

setMethod("subsetRows", "ClusteredScoresList",
    function(x, i = NULL)
{
    if(is.null(i))
        stop("No row indices given to subset by.")

    old.ranges <- lapply(x@scores, range, na.rm = TRUE)
    new("ClusteredScoresList", names = x@names,
        scores = lapply(x@scores, function(y) y[i, ]),
	anno = x@anno[i], up = x@up, down = x@down, dist = x@dist,
	freq = x@freq, s.width = x@s.width, cluster.id = x@cluster.id[i],
	expr = x@expr[i], expr.name = x@expr.name, sort.data = x@sort.data[i],
        sort.name = x@sort.name, .old.ranges = old.ranges, .samp.info = x@.samp.info)
})

setGeneric("clusters", function(x, ...)
           {standardGeneric("clusters")})
setMethod("clusters", "ClusteredScoresList",
    function(x)
{
    x@cluster.id
})

setClass("GCAdjustParams", representation(
                                    genome = "BSgenome",
                                    mappability = "BSgenome",
                                    min.mappability = "numeric",
				    n.bins = "numeric",
                                    min.bin.size = "numeric",
				    poly.degree = "numeric",
                                    ploidy = "numeric")
)

# Constructor
setGeneric("GCAdjustParams", function(genome, mappability, ...)
           {standardGeneric("GCAdjustParams")})
setMethod("GCAdjustParams", c("BSgenome", "BSgenome"),
    function(genome, mappability, min.mappability, n.bins = NULL, min.bin.size = 1,
             poly.degree = NULL, ploidy = 1)
{
    if(is.null(min.mappability))
        stop("Minimum mappability of counting windows to keep not given.")

    if(is.null(n.bins))
        stop("Number of GC bins to bin counts into not given.")

    if(is.null(poly.degree))
        stop("Polynomial degree not given.")

    new("GCAdjustParams", genome = genome, mappability = mappability, min.mappability = min.mappability,
        n.bins = n.bins, min.bin.size = min.bin.size, poly.degree = poly.degree,
        ploidy = ploidy)
})

# container for output of regionStats()    
setClass("RegionStats", representation("list"))

setMethod("show", "RegionStats", function(object) {
  cat("Object of class 'RegionStats'.\n")
  cat("Results for: ", paste(names(object$regions),collapse=" "), "\n")
  cat("Names:", paste(names(object),collapse=" "), "\n")
})

# container for output of ChromaBlocks()
setClass("ChromaResults",
    representation(
        blocks="GRanges", 
        regions="RangesList",
        FDRTable="matrix",
        cutoff="numeric"
    )
)

setMethod("show", "ChromaResults", function(object) {
  cat("Object of class 'ChromaResults'.\n")
  cat(sum(sapply(object@regions, length)), "regions found with using a cutoff of", object@cutoff, "\n")
})

#ChromaResults Generics
setGeneric("blocks", function(x) standardGeneric("blocks"))
setGeneric("regions", function(x) standardGeneric("regions"))
setGeneric("FDRTable", function(x) standardGeneric("FDRTable"))
setGeneric("cutoff", function(x) standardGeneric("cutoff"))

#ChromaResults Accessors
setMethod("blocks", "ChromaResults", function(x) x@blocks)
setMethod("regions", "ChromaResults", function(x) x@regions)
setMethod("FDRTable", "ChromaResults", function(x) x@FDRTable)
setMethod("cutoff", "ChromaResults", function(x) x@cutoff)

