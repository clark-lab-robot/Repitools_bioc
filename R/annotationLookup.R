setGeneric("annotationBlocksLookup", function(x, anno, ...)
           {standardGeneric("annotationBlocksLookup")})
setGeneric("annotationLookup", function(x, anno, ...)
           {standardGeneric("annotationLookup")})
setGeneric("annotationBlocksCounts", function(x, anno, ...)
           {standardGeneric("annotationBlocksCounts")})
setGeneric("annotationCounts", function(x, anno, ...)
           {standardGeneric("annotationCounts")})

setMethod("annotationBlocksLookup", c("data.frame", "GRanges"),
    function(x, anno, verbose = TRUE)
{
    if("index" %in% colnames(x)) p.inds <- x$index else p.inds <- 1:nrow(x) 
    probesGR <- GRanges(x$chr, IRanges(x$position, width = 1))
    g.names <- .getNames(anno)

    if(verbose) message("Processing mapping between probes and features.")
    mapping <- suppressWarnings(findOverlaps(anno, probesGR)@matchMatrix)
    inds <- split(p.inds[mapping[, 2]], factor(mapping[, 1],
                                levels = 1:length(g.names)))
    pos.list <- split(start(probesGR[mapping[, 2]]), factor(mapping[, 1],
                                             levels = 1:length(g.names)))
    
    offs <- mapply(function(u, v) u - v, pos.list, start(anno),
                   SIMPLIFY = FALSE)
    offs <- mapply(function(u, v)
                   {
                       names(v) <- u
                       v
                   }, inds, offs, SIMPLIFY = FALSE)
    names(inds) <- names(offs) <- g.names
    if(verbose) message("Mapping done.")

    list(indexes = inds, offsets = offs)
})

setMethod("annotationBlocksLookup", c("data.frame", "data.frame"),
    function(x, anno, ...)
{
    col.missing <- setdiff(c("chr", "start", "end"), colnames(anno))
    if(length(col.missing) > 0)
	stop("Columns ", paste(col.missing, collapse = ", "),
             " of annotation are not present.")

    annotationBlocksLookup(x, annoDF2GR(anno), ...)
})

setMethod("annotationLookup", c("data.frame", "GRanges"),
    function(x, anno, up, down, ...)
{
    invisible(.validate(anno, up, down))

    blocksGR <- featureBlocks(anno, up, down)
    annot <- annotationBlocksLookup(x, blocksGR, ...)

    pos <- as.logical(strand(anno) == '+')
    annot$offsets[pos] <- lapply(annot$offsets[pos], function(z) z - up)
    annot$offsets[!pos] <- lapply(annot$offsets[!pos], function(z) rev(down - z))
    annot$indexes[!pos] <- lapply(annot$indexes[!pos], rev)
    names(annot$offsets) <- names(annot$indexes) <- .getNames(anno)

    annot
})

setMethod("annotationLookup", c("data.frame", "data.frame"),
    function(x, anno, ...)
{
    col.missing <- setdiff(c("chr", "start", "end", "strand"), colnames(anno))
    if(length(col.missing) > 0)
	stop("Columns ", paste(col.missing, collapse = ", "),
             " of annotation are not present.")

    annotationLookup(x, annoDF2GR(anno), ...)
})

setMethod("annotationBlocksCounts", c("GRanges", "GRanges"),
    function(x, anno, seq.len = NULL, verbose = TRUE)
{
    require(GenomicRanges)

    f.names <- .getNames(anno)
    if(!is.null(seq.len))
    {
        if(verbose)
            message("Extending all reads to fragment length.")
        x <- resize(x, seq.len)
        if(verbose)
            message("Read extension complete.\nCounting started.")
    }
    counts <- matrix(countOverlaps(anno, x))
    rownames(counts) <- f.names
    if(verbose)
    	message("Counting successful.")
    
    counts
})

setMethod("annotationBlocksCounts", c("GRangesList", "GRanges"),
    function(x, anno, ...)
{
    require(GenomicRanges)

    f.names <- .getNames(anno)
    counts <- IRanges::lapply(x, function(z)
                  annotationBlocksCounts(z, anno, ...))
    counts <- do.call(cbind, counts)
    if(!is.null(names(x))) colnames(counts) <- names(x)
    rownames(counts) <- f.names

    counts
})

setMethod("annotationBlocksCounts", c("character", "GRanges"),
    function(x, anno, ...)
{
    f.names <- .getNames(anno)
    counts <- lapply(x, function(z) annotationBlocksCounts(BAM2GRanges(z), anno, ...))
    counts <- do.call(cbind, counts)
    if(!is.null(names(x))) colnames(counts) <- names(x)
    rownames(counts) <- f.names

    counts
})

setMethod("annotationBlocksCounts", c("ANY", "data.frame"),
    function(x, anno, ...)
{
    annotationBlocksCounts(x, annoDF2GR(anno), ...)
})

setMethod("annotationCounts", c("ANY", "GRanges"),
    function(x, anno, up, down, ...)
{
    invisible(.validate(anno, up, down))
    blocksGR <- featureBlocks(anno, up, down)

    annotationBlocksCounts(x, blocksGR, ...)
})

setMethod("annotationCounts", c("ANY", "data.frame"),
    function(x, anno, ...)
{
    annotationCounts(x, annoDF2GR(anno), ...)
})
