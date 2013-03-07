# A collection of variables that describe where the sampling will happen.
# An S4 class, so that they are created once, then dispatched on when
# required for multiple samples.

setClass(".CoverageSamples",
         representation(
                        pos.labels = "ANY", # character or numeric.
                        cvg.samps = "GRanges",
                        marks.samps.map = "ANY" # NULL or BSgenome.
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
                                    cluster.id = "factor",
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
    function(x, anno = x@anno, scores = tables(x), expr = NULL, expr.name = NULL, cluster.id, sort.data = NULL,
             sort.name = NULL)
{
	new("ClusteredScoresList", names = x@names, scores = scores, anno = anno,
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
    function(genome, mappability, min.mappability, n.bins = NULL, min.bin.size = 2,
             poly.degree = NULL, ploidy = 1)
{
    if(is.null(min.mappability))
        stop("Minimum mappability of counting windows to keep not given.")

    if(is.null(n.bins))
        stop("Number of GC bins to bin counts into not given.")

    new("GCAdjustParams", genome = genome, mappability = mappability, min.mappability = min.mappability,
        n.bins = n.bins, min.bin.size = min.bin.size, poly.degree = poly.degree,
        ploidy = ploidy)
})

setClass("CopyEstimate", representation(
                                    windows = "GRanges",
                                    unadj.CN = "matrix",
                                    unadj.CN.seg = "GRangesList",
                                    type = "character")
)
setGeneric("CopyEstimate", function(windows, unadj.CN, unadj.CN.seg, type)
           {standardGeneric("CopyEstimate")})

setMethod("CopyEstimate", c("GRanges", "matrix", "GRangesList", "character"),
    function(windows, unadj.CN, unadj.CN.seg, type = "relative")
{
    new("CopyEstimate", windows = windows, unadj.CN = unadj.CN, unadj.CN.seg = unadj.CN.seg, type = type)
})

setMethod("show", "CopyEstimate", function(object) {
    cat("Object of class 'CopyEstimate'.\n")
    cat("Windows:\n")
    print(object@windows)
    cat("Unadjusted Copy Number (first 6):\n")
    print(head(object@unadj.CN))
    cat("Data Type: ", object@type, "\n", sep = '')
    if(length(object@unadj.CN.seg) > 0)
    {
        cat("Segmented Copy Number Estimates:\n")
        print(object@unadj.CN.seg)
    }
})

setClass("AdjustedCopyEstimate", representation(
                                    ploidy = "numeric",
                                    models = "list",
                                    adj.CN = "matrix",
                                    adj.CN.seg = "GRangesList"),
                                contains = "CopyEstimate"
)
setGeneric("AdjustedCopyEstimate", function(ploidy, windows, mappability, gc, unadj.CN, models, adj.CN, type)
           {standardGeneric("AdjustedCopyEstimate")})

setMethod("AdjustedCopyEstimate", c("numeric", "GRanges", "numeric", "numeric", "matrix", "list", "matrix", "character"),
    function(ploidy, windows, mappability, gc, unadj.CN, models, adj.CN, type)
{
    if(length(ploidy) < ncol(unadj.CN))
        ploidy <- rep(ploidy, ncol(unadj.CN))

    elementMetadata(windows) <- DataFrame(mappability, GC = gc)

    new("AdjustedCopyEstimate", ploidy = ploidy, windows = windows, unadj.CN = unadj.CN,
        models = models, adj.CN = adj.CN, type = type)
})

setMethod("show", "AdjustedCopyEstimate", function(object) {
    cat("Object of class 'AdjustedCopyEstimate'.\n")
    cat("Genome Ploidy: ", paste(object@ploidy, collapse = ", "), '\n', sep = '')
    cat("Windows:\n")
    print(object@windows)
    cat("Unadjusted Copy Number (first 6):\n")
    print(head(object@unadj.CN))
    if(length(object@unadj.CN.seg) > 0)
    {
        cat("Segmented Unadjusted Copy Number Estimates:\n")
        print(object@unadj.CN.seg)
    }
    cat("Model Fits:\n")
    print(object@models)
    cat("Adjusted Copy Number (first 6):\n")
    print(head(round(object@adj.CN, 2)))
    if(length(object@adj.CN.seg) > 0)
    {
        cat("Segmented Adjusted Copy Number Estimates:\n")
        print(object@adj.CN.seg)
    }
    cat("Data Type: ", object@type, "\n", sep = '')
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


# container for abcdDNA stuff
setClass("QdnaData",representation("list"))

#    counts="matrix", 
#    regions="GRanges",
#    offsets="matrix",
#    neutral="logical",
#    design="matrix",
#    cnv="matrix"
#  )
#)

#setMethod("show", "QdnaData", function(object) {
#	cat("Object of class 'QdnaData'.\n")
#	cat(nrow(object@counts), "regions for", ncol(object@counts),"samples\n")
#	cat("Slots:", slotNames(object),"\n")
#})

setMethod("show", "QdnaData", function(object) {
	cat("Object of class 'QdnaData'.\n")
	cat(nrow(object$DGEList$counts), "regions for", ncol(object$DGEList$counts),"samples\n")
	cat("Slots:", names(object),"\n")
})


#cat("\nCounts:\n")
#	print(head(object@counts))
#	cat("\nRegions:\n")
#	print(object@regions)
#	cat("\nOffsets:\n")
#	print(head(object@offsets))

#setGeneric("QdnaData", function(counts,regions,design,cnv,offsets,neutral) { standardGeneric("QdnaData") })
#setMethod("QdnaData", c("GRanges", "matrix", "GRangesList", "character"),
#function(counts,regions,design,cnv,offsets,neutral) {
#    new("CopyEstimate", windows = windows, unadj.CN = unadj.CN, unadj.CN.seg = unadj.CN.seg, type = type)
#})


######################################################################
## Definition of the "BayMethList" class
######################################################################

setClass("BayMethList", representation(
    windows="GRanges", ## Info to which genomic windows the data belong to
    control="matrix",  ## SssI control 
    sampleInterest="matrix", ## Sample of interest 
    cpgDens="numeric",  ## CpG density
    f="matrix", ## Normalizing offset (possibly including CN-variations)
    priorTab="list", ## List of prior parameters for each sample
    methEst="list", ## List to save the results
    maskEmpBayes="logical" ## indicate which bins should be masked out in the empirical Bayes
    ))


######################################################################
## Constructor for the "BayMethList" class
######################################################################
BayMethList <- function(windows, control, sampleInterest, cpgDens,
f=matrix(), priorTab=list(), methEst=list(), maskEmpBayes=logical())
{
    if((class(windows) != "GRanges") || (class(control) != "matrix") ||
        (class(sampleInterest) != "matrix") || 
        (class(cpgDens) != "numeric") ){
            stop("\n\n\t `windows' must be of class `GRanges', `control' and
                `sampleInterest' of class `matrix' and `cpgDens' a `numeric'
                vector\n\n")
    }
    if(class(f) != "matrix"){
        stop("\n\n\t 'f' must be a matrix\n\n")
    }
    if(class(priorTab) != "list"){
        stop("\n\n\t 'priorTab' must be a list\n\n")
    }
    if(class(methEst) != "list"){
        stop("\n\n\t 'methEst' must be a list\n\n")
    }
    if(class(maskEmpBayes) != "logical"){
        stop("\n\n\t 'maskEmpBayes' must be a logical vector\n\n")
    }
    ## get the object dimensions
    na <- length(windows)
    nc <- nrow(control)
    ns <- nrow(sampleInterest)
    ncp <- length(cpgDens)
    nf <- nrow(f)
    nm <- length(maskEmpBayes)

    ## either control is a matrix with one column or 
    ## with the same number of columns as sampleInterest
    if((ncol(control) != 1) && (ncol(control) != ncol(sampleInterest))){
        stop("\n\n\tNumber of SssI-controls is not correct. 
                    Either only one control or as many as sample of 
                    interests need to be provided.\n\n")
    }
    if(length(unique(c(na, nc, ns, ncp))) != 1){
        stop("\n\n\tThe annotation matrix, SssI control matrix, 
                    sample of interest matrix and CpG density have not 
                    the same length.\n\n")
    }
    if((nf != 1)  && (nf != na)){
        stop("\n\n\tThe number of offsets per sample must be either one or be
equal to the length of the annotation matrix.\n\n")
    }
    if((nm != 0)  && (nm != na)){
        stop("\n\n\tThe logical vector to mask bins out from the empirical Bayes approach must be of the same lenth as the number of bins.\n\n")
    }
    if(nm == 0){
        maskEmpBayes <- rep(FALSE, na)
    }
    new("BayMethList", windows=windows, control=control, 
        sampleInterest=sampleInterest, cpgDens=cpgDens, f=f, 
        priorTab=priorTab, methEst=methEst, maskEmpBayes=maskEmpBayes)
}


######################################################################
## Show function for the "BayMethList" class
######################################################################

setMethod("show", "BayMethList", function(object) {
    cat("Object of class 'BayMethList'.\n\n")
    cat("- Genomic windows have width:", unique(width(object@windows)), "\n\n")
    print(seqnames(object@windows))
    cat("\n- Number of control samples:", ncol(object@control), 
        "\n")
    #print(apply(object@control, 2, summary))
    cat("\n- Number of samples of interest:",
        ncol(object@sampleInterest), "\n")
    #print(apply(object@sampleInterest, 2, summary))
    #cat("\nSummary information for the CpG density:\n")
    #print(summary(object@cpgDens))
    if((dim(object@f) != c(1,1)) || (!is.na(object@f[1,1]))){
        cat("\n- Slot for normalizing offset is filled.\n")
    } else {
        cat("\n- Slot for normalizing offset is empty.\n") 
    }
    if(length(object@priorTab) > 0){
        cat("\n- Prior parameters are available.\n")
    } else {
        cat("\n- Prior parameters are NOT available.\n")        
    }
    if(length(object@methEst) > 0){
        cat("\n- Methylation estimates are available.\n")
    } else {
        cat("\n- Methylation estimates are NOT available.\n")
    }
    if(length(object@maskEmpBayes) > 0){
        cat("\n- ", sum(object@maskEmpBayes), " bins are masked out.\n")
    }
})

######################################################################
## Access functions for the "BayMethList" class
######################################################################

if(!isGeneric("[")) setGeneric("[", function(object) standardGeneric("["))
setMethod("[", "BayMethList",
    function(x, i) {

    message("\n\n\tCAUTION: Slots 'f', 'priorTab' and'methEst' ", 
        "do not change when taking the subset!\n\n")

    BayMethList(windows=x@windows[i], 
        control=matrix(x@control[i,], ncol=ncol(x@control)),
        sampleInterest=matrix(x@sampleInterest[i,], 
            ncol=ncol(x@sampleInterest)), 
        cpgDens=x@cpgDens[i],
        f=x@f,
        priorTab=x@priorTab,
        methEst=x@methEst,
        maskEmpBayes=x@maskEmpBayes[i])
})
setMethod("length", "BayMethList",
    function(x) {
        length(x@windows)
})
if(!isGeneric("windows")) setGeneric("windows", 
    function(object) standardGeneric("windows"))
setMethod("windows", "BayMethList", 
    function(object) {
        object@windows
})
if(!isGeneric("control")) setGeneric("control", 
    function(object) standardGeneric("control"))
setMethod("control", "BayMethList", 
    function(object) {
        object@control
})  
if(!isGeneric("sampleInterest")) setGeneric("sampleInterest", 
    function(object) standardGeneric("sampleInterest"))
setMethod("sampleInterest", "BayMethList", 
    function(object) {
        object@sampleInterest
})
if(!isGeneric("cpgDens")) setGeneric("cpgDens", 
    function(object) standardGeneric("cpgDens"))
setMethod("cpgDens", "BayMethList", 
    function(object) {
        object@cpgDens
})
if(!isGeneric("fOffset")) setGeneric("fOffset", 
    function(object) standardGeneric("fOffset"))
setMethod("fOffset", "BayMethList", 
    function(object) {
        object@f
})
if(!isGeneric("priorTab")) setGeneric("priorTab", 
    function(object) standardGeneric("priorTab"))
setMethod("priorTab", "BayMethList", 
    function(object) {
        object@priorTab
})
if(!isGeneric("methEst")) setGeneric("methEst", 
    function(object) standardGeneric("methEst"))
setMethod("methEst", "BayMethList", 
    function(object) {
        object@methEst
})
if(!isGeneric("maskEmpBayes")) setGeneric("maskEmpBayes", 
    function(object) standardGeneric("maskEmpBayes"))
setMethod("maskEmpBayes", "BayMethList", 
    function(object) {
        object@maskEmpBayes
})


## Determine some lengths
if(!isGeneric("ncontrol")) setGeneric("ncontrol", 
    function(x) standardGeneric("ncontrol"))
setMethod("ncontrol", "BayMethList", function(x) ncol(x@control))
if(!isGeneric("nsampleInterest")) setGeneric("nsampleInterest", 
    function(x) standardGeneric("nsampleInterest"))
setMethod("nsampleInterest", "BayMethList", function(x) ncol(x@sampleInterest))

######################################################################
## Replace functions for the "BayMethList" class
######################################################################

setGeneric("windows<-", function(x, value) standardGeneric("windows<-"))
setReplaceMethod("windows", "BayMethList", function(x, value) {
    x@windows <- value
    x
})
setGeneric("control<-", function(x, value) standardGeneric("control<-"))
setReplaceMethod("control", "BayMethList", function(x, value) {
    x@control <- value
    x
})
setGeneric("sampleInterest<-", 
function(x, value) standardGeneric("sampleInterest<-"))
setReplaceMethod("sampleInterest", "BayMethList", function(x, value) {
    x@sampleInterest <- value
    x
})
setGeneric("cpgDens<-", function(x, value) standardGeneric("cpgDens<-"))
setReplaceMethod("cpgDens", "BayMethList", function(x, value) {
    x@cpgDens <- value
    x
})
setGeneric("fOffset<-", function(x, value) standardGeneric("fOffset<-"))
setReplaceMethod("fOffset", "BayMethList", function(x, value) {
    
    if(class(value) != "matrix"){
        stop("The offset must be of class matrix with the same number of
columns as sample of interests.")
    }
   
    if((nrow(value) != 1) && (nrow(value) != length(x))){
        stop("\n\n\tThe number of offsets per sample must be either one or be
equal to the length of the annotation matrix.\n\n")
    }

    x@f <- value
    x
})
setGeneric("priorTab<-", function(x, value) standardGeneric("priorTab<-"))
setReplaceMethod("priorTab", "BayMethList", function(x, value) {
    x@priorTab <- value
    x
})
setGeneric("methEst<-", function(x, value) standardGeneric("methEst<-"))
setReplaceMethod("methEst", "BayMethList", function(x, value) {
    x@methEst <- value
    x
})
setGeneric("maskEmpBayes<-", function(x, value) standardGeneric("maskEmpBayes<-"))
setReplaceMethod("maskEmpBayes", "BayMethList", function(x, value) {
    x@maskEmpBayes <- value
    x
})




