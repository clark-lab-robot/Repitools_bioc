setGeneric("featureScores", signature = c("x", "anno"), function(x, anno, ...)
                                           {standardGeneric("featureScores")})
setGeneric(".featureScores", signature = c("x", "y"), function(x, y, ...)
                                           {standardGeneric(".featureScores")})

setClassUnion(".SequencingData", c("character", "GRanges", "GRangesList"))

setMethod(".featureScores", c("GRanges", ".CoverageSamples"),
    function(x, y, anno, up, down, dist, freq, s.width, use.strand = FALSE, verbose)
{
    require(GenomicRanges)

    # Unpack variables in y.
    pos.labels <- y@pos.labels
    cvg.samps <- y@cvg.samps
    if(use.strand == FALSE)
        strand(cvg.samps) <- '*'

    # Only use sequencing data on annotated chromosomes.
    x <- x[seqnames(x) %in% seqlevels(anno)]
    seqlevels(x) <- seqlevels(anno)

    # Qualitatively near identical to running mean smoothing.
    if(verbose) message("Extending all reads to smoothing width.")
    seqlengths(x) <- rep(NA, length(seqlengths(x)))
    if(!is.null(s.width))
        x <- resize(x, s.width)

    # Get coverage.
    if(verbose) message("Calculating coverage at sample points.")
    # Scale coverages for total reads.
    cvg.mat <- matrix(countOverlaps(cvg.samps, x) / length(x),
                      ncol = length(pos.labels),
                      byrow = TRUE)

    # Precision sometimes means 0 is represented as very small negative numbers.
    cvg.mat[cvg.mat < 0] = 0
    	
    colnames(cvg.mat) <- pos.labels
    rownames(cvg.mat) <- .getNames(anno)

    new("ScoresList", names = "Undefined", scores = list(cvg.mat), anno = anno,
         up = up, down = down, dist = dist, freq = freq, s.width = s.width,
         .samp.info = y)
})

setMethod(".featureScores", c("GRangesList", ".CoverageSamples"),
    function(x, y, anno, up, down, dist, freq, s.width, ..., verbose)
{
    if(length(s.width) == 1)
        s.width <- rep(s.width, length(x))
    scores <- mapply(function(z, i)
	           {
                        if(verbose && !is.null(names(x)))
                            message("Processing sample ", names(x)[i])
		   	.featureScores(z, y, anno, up, down, dist,
                                        freq, s.width[i], ..., verbose = verbose)
		   }, x, IntegerList(as.list(1:length(x))), SIMPLIFY = FALSE)

    if(!is.null(names(x)))
	names <- names(x)
    else
	names <- unname(sapply(scores, names))
    new("ScoresList", names = names, anno = anno, scores = unname(sapply(scores, tables)),
                            up = up, down = down, dist = dist, freq = freq,
                            s.width = s.width, .samp.info = y)
})

setMethod(".featureScores", c("character", ".CoverageSamples"),
    function(x, y, anno, up, down, dist, freq, s.width, ..., verbose)
{
    if(length(s.width) == 1)
        s.width <- rep(s.width, length(x))

    scores <- mapply(function(z, i)
	           {
                        if(verbose && !is.null(names(x)))
                            message("Processing sample ", names(x)[i])
		   	
		   	.featureScores(BAM2GRanges(z), y, anno, up, down, dist, freq,
                                        s.width[i], ..., verbose = verbose)
		   }, x, 1:length(x), SIMPLIFY = FALSE)
    if(!is.null(names(x)))
	names <- x
    else
	names <- unname(sapply(scores, names))
    new("ScoresList", names = names, anno = anno, scores = unname(sapply(scores, tables)),
                            up = up, down = down, dist = dist, freq = freq,
                            s.width = s.width, .samp.info = y)
})

setMethod(".featureScores", c(".SequencingData", "GRanges"),
    function(x, y, up, down, dist = c("base", "percent"), freq, s.width, ...,
             verbose = TRUE)
{
    require(IRanges)
    dist <- match.arg(dist)

    str <- strand(y)
    st <- start(y)
    en <- end(y)	
    wd <- width(y) 
    pos <- str == '+'

    if(verbose) message("Calculating sampling positions.")
    cov.winds <- featureBlocks(y, up, down, dist, keep.strand = TRUE)

    posns <- seq(-up, down, freq)
    if(dist == "percent")
	pos.labels <- paste(posns, '%')
    else
	pos.labels <- posns
    n.pos <- length(posns)

    # Make ranges for each sample point.
    cvg.samps <- rep(cov.winds, each = n.pos)
    if(dist == "percent")
        gap.size <- width(cvg.samps) / (n.pos - 1)
    else
        gap.size <- rep(freq, length(cvg.samps))

    ranges(cvg.samps) <- IRanges(start = as.numeric(
                              ifelse(strand(cvg.samps) %in% c('+', '*'), 
                                     start(cvg.samps) + 0:(n.pos - 1) * gap.size,
                                     end(cvg.samps) - 0:(n.pos - 1) * gap.size)
                                                   ),
                                     width = 1)

    samp.info <- new(".CoverageSamples", pos.labels = pos.labels, cvg.samps = cvg.samps)

    .featureScores(x, samp.info, y, up, down, dist, freq, s.width, ..., verbose = verbose)
})

setMethod(".featureScores", c("matrix", "GRanges"),
    function(x, y, up, down, p.anno, mapping = NULL, freq, log2.adj = TRUE, verbose = TRUE)
{
    if(is.null(p.anno))
        stop("Probe annotation not given.")
    if(is.null(freq))
        stop("Sampling frequency not given.")

    if(is.null(mapping))
    {
        if("index" %in% colnames(x)) p.inds <- x$index else p.inds <- 1:nrow(x)    
        ind.col <- colnames(p.anno) == "index"
        mapping <- annotationLookup(p.anno[, !ind.col], y, up, down, verbose)
        p.used <- unique(unlist(mapping$indexes, use.names = FALSE))
        p.anno <- p.anno[p.used, ]
        mapping <- annotationLookup(p.anno[, !ind.col], y, up, down, verbose)
        intens <- x[p.anno$index, ]
    } else {
        intens <- x
    }

    if(log2.adj) intens <- log2(intens)

    points.probes <- makeWindowLookupTable(mapping$indexes, mapping$offsets,
                     starts = seq(-up, down - freq, freq),
                     ends = seq(-up + freq, down, freq))

    points.intens <- lapply(1:ncol(x), function(z)
                     {
                         .scoreIntensity(points.probes, intens[, z], returnMatrix = TRUE)
                     })    

    new("ScoresList", names = colnames(x), anno = y, scores = points.intens,
                            up = up, down = down, dist = "base", freq = freq,
                            s.width = NULL)
})

setMethod(".featureScores", c("AffymetrixCelSet", "GRanges"),
    function(x, y, p.anno = NULL, mapping = NULL, chrs = NULL, ...)
{
    require(aroma.affymetrix)

    if(is.null(mapping) && is.null(p.anno))
        p.anno <- getProbePositionsDf(getCdf(x), chrs, verbose = verbose)

    intens <- extractMatrix(x, cells = p.anno$index, verbose = verbose)
    p.anno$index <- 1:nrow(p.anno)

    .featureScores(intens, y, p.anno = p.anno, mapping = mapping, ...)
})

setMethod("featureScores", c("ANY", "GRanges"), function(x, anno,
           up = NULL, down = NULL, ...)
{
    invisible(.validate(anno, up, down))
    .featureScores(x, anno, up = up, down = down, ...)
})

setMethod("featureScores", c("ANY", "data.frame"),
    function(x, anno, ...)
{
    featureScores(x, annoDF2GR(anno), ...)
})
