setGeneric("blocksStats", function(x, anno, ...){standardGeneric("blocksStats")})
setGeneric(".blocksStats", function(x, anno, ...){standardGeneric(".blocksStats")})

setMethod(".blocksStats", c("matrix", "GRanges"),
    function(x, anno, up, down, p.anno = NULL, mapping = NULL, log2.adj = TRUE,
             design, robust = 0, p.adj = "fdr", verbose = TRUE)
{
    if(nrow(design) != ncol(x))
        stop("The number of rows in the design matrix does not equal
              the number of columns in the probes data matrix.\n")

    if(nrow(x) != nrow(p.anno))
        stop("Number of rows in probe annotation and intensity matrix differ.\n")
		
    used <- which(rowSums(design != 0) > 0)
	
    if(is.null(mapping))
    {	
        ind.col <- colnames(p.anno) == "index"
        if(all(ind.col == FALSE))
            stop("'index' column missing from probe annotation.\n")
        if(!"name" %in% colnames(elementMetadata(anno)))
            stop("'name' column missing from annotation.\n")

        # Run lookup twice. First to get a list of smaller list of probes to use.
	if(is.null(up))
            mapping <- annotationBlocksLookup(p.anno[, !ind.col], anno,
                                              verbose = verbose)
        else
            mapping <- annotationLookup(p.anno[, !ind.col], anno, up, down,
                                        verbose = verbose)
	p.used <- unique(unlist(mapping$indexes, use.names = FALSE))
	p.anno <- p.anno[p.used, ]
	if(is.null(up))
            mapping <- annotationBlocksLookup(p.anno[, !ind.col], anno,
                                              verbose = verbose)
        else
            mapping <- annotationLookup(p.anno[, !ind.col], anno, up, down,
                                        verbose = verbose)
    }

    if(log2.adj == TRUE)
        intens <- log2(x[p.anno$index, used])
    else
        intens <- x[p.anno$index, used]

    diffs <- intens %*% design[used, , drop = FALSE]

    means <- t.stats <- matrix(NA, nrow = length(anno), ncol = ncol(diffs),
                              dimnames = list(NULL, colnames(design)))
    df <- rep(0, length(anno))

    for(i in 1:length(anno))
    {
        pr.inds <- mapping[[1]][[i]]
        pr.inds <- unique(pr.inds[!is.na(pr.inds)])
        if(length(pr.inds) < 2)
            next
        for(j in 1:ncol(diffs))
        {
            if(robust && length(pr.inds) >= robust)
            {
                require(MASS)
                r.model <- summary(rlm(diffs[pr.inds, j] ~ 1))
                t.stats[i, j] <- r.model$coef[3]
                means[i, j] <- r.model$coef[1]
                df[i] <- r.model$df[2]
            } else {
                    t.vals <- t.test(diffs[pr.inds, j])
                    t.stats[i, j] <- t.vals$statistic
                    means[i, j] <- t.vals$estimate
                    df[i] <- t.vals$parameter
            }
        }
    }

    p.vals <- 2 * pt(-abs(t.stats), df)
  
    # Adjust p-values for multiple testing.
    adj.p <- p.vals
    for(i in 1:ncol(adj.p))
        adj.p[, i] <- p.adjust(p.vals[, i], method = p.adj)
  
    results <- data.frame(df = df, mean.diff = means, t.stats = t.stats,
                          p.vals = p.vals, adj.p.vals = adj.p)

    if(is.null(up))
        colnames(results)[1] <- "df"        
    else
        colnames(results)[1] <- paste("df", up, down, sep = '.')
    
    if(ncol(design) == 1)    
        colnames(results)[2:5] <- paste(c("meandiff", "t.stats", "p.vals",
                                  "adj.p.vals"), gsub(".[1-9]$", '',
                                  colnames(results)[2:5]), sep = '.')

    cbind(annoGR2DF(anno), results)
})

setMethod(".blocksStats", c("AffymetrixCelSet", "GRanges"),
    function(x, anno, up, down, p.anno = NULL, chrs = NULL, mapping = NULL,
             design = NULL, ...)
{
    if(is.null(design))
        stop("No design matrix given.")
    require(aroma.affymetrix)

    if(nrow(design) != nbrOfArrays(x))
        stop("The number of rows in the design matrix does not equal the number
              of arrays in the CEL set.\n")

    used <- which(rowSums(design != 0) > 0)
    x <- extract(x, used, verbose = verbose)

    if(is.null(mapping) && is.null(p.anno))
        p.anno <- getProbePositionsDf(getCdf(x), chrs, verbose = verbose)
    intens <- extractMatrix(x, cells = p.anno$index, verbose = verbose)
    p.anno$index <- 1:nrow(p.anno)

    .blocksStats(intens, anno, up, down, p.anno, mapping,
                 design = design[used, , drop = FALSE], ...)
})

setMethod(".blocksStats", c("GRangesList", "GRanges"),
    function(x, anno, up, down, seq.len = NULL, design = NULL, lib.size = "lane",
             Acutoff = NULL, p.adj = "fdr", verbose = TRUE)
{
    if(is.null(design))
        stop("No design matrix given.")

    require(edgeR)

    if(lib.size == "ref" && is.null(Acutoff))
        stop("Must give value of Acutoff if using \"ref\" normalisation.\n")

    if (is.null(up))
        counts <- annotationBlocksCounts(x, anno, seq.len, verbose)
    else
        counts <- annotationCounts(x, anno, up, down, seq.len, verbose)
	
    lib.sizes <- switch(lib.size,
                        lane = elementLengths(x),
                        blocks = colSums(counts),
                        ref = colSums(counts) * calcNormFactors(counts,
                              Acutoff = Acutoff))
    annoDF <- annoGR2DF(anno)
    results <- cbind(annoDF, counts)
    for (i in 1:ncol(design))
    {
	if (verbose == TRUE)
            message("Processing column ", i, " of design matrix.")
	stopifnot(sum(design[,i] ==  1) > 0,
                  sum(design[,i] == -1) > 0,
                  all(design[,i] %in% c(-1, 0, 1)))
	used <- design[, i] != 0
	dge <- DGEList(counts = counts[,used],
                     group = as.character(design[used, i]),
                     lib.size = lib.sizes[used])
	dge <- estimateCommonDisp(dge)
	pseudo.counts <- dge$pseudo.alt
	colnames(pseudo.counts) <- paste(colnames(pseudo.counts), "pseudo",
                                         sep = '_')
	results <- cbind(results, pseudo.counts)
	de.stats <- exactTest(dge, pair = c("-1", "1"))$table
        adj.p.vals <- p.adjust(de.stats[, 3], p.adj)
        de.stats <- cbind(de.stats, adj.p.vals)
	colnames(de.stats) <- paste(colnames(de.stats), colnames(design)[i],
                                    sep = "_")

	results <- cbind(results, de.stats)
    }
    results
})

setMethod(".blocksStats", c("character", "GRanges"),
    function(x, anno, ...)
{
    .blocksStats(BAM2GRangesList(x), anno, ...)
})

setMethod("blocksStats", c("ANY", "GRanges"),
    function(x, anno, up = NULL, down = NULL, ...)
{
    if(!is.null(up) && !is.null(down))
        invisible(.validate(anno, up, down))
    .blocksStats(x, anno, up, down, ...)
})

setMethod("blocksStats", c("ANY", "data.frame"),
    function(x, anno, ...)
{
    blocksStats(x, annoDF2GR(anno), ...)
})
