setGeneric("closeClusterFeatures", function(c.list, ...){standardGeneric("closeClusterFeatures")})

.clusterBlocksInds <- function(cluster.IDs.chr, n.consec)
{
    id.levels <- sort(unique(cluster.IDs.chr))
    blocks.inds.chr <- lapply(id.levels, function(x)
    {
        which.id <- IRanges(which(cluster.IDs.chr == x), width = 1)
        which.in.blocks <- reduce(which.id)
        which.in.blocks <- which.in.blocks[width(which.in.blocks) >= 3]
        which.dispersed <- reduce(which.id, min.gapwidth = 2)
        unique(which.dispersed[subjectHits(findOverlaps(which.in.blocks, which.dispersed))])
    })
    names(blocks.inds.chr) <- id.levels
    blocks.inds.chr
}

.clusterBlocks <- function(tx.df, cluster.IDs, name.col, ...)
{
    anno.cols <- colnames(tx.df) %in% c("gene", "gene.symbol")
    chr.blocks.inds <- tapply(cluster.IDs, tx.df[, "seqnames"], function(x) .clusterBlocksInds(x, ...))
    blocks.inds <- unlist(mapply(function(clust.blocks, chr)
    {
        which.on.chr <- which(tx.df[, "seqnames"] == chr)
        lapply(clust.blocks, function(blocks)
        {
            sapply(as.list(blocks), function(block)
                                    {
                                        feat.inds <- which.on.chr[block]
                                        feat.locs <- GRanges(tx.df[feat.inds, "seqnames"],
                                                             IRanges(tx.df[feat.inds, "start"], tx.df[feat.inds, "end"]), 
                                                             tx.df[feat.inds, "strand"])
                                        values(feat.locs) <- DataFrame(feature.id = tx.df[feat.inds, 5 + name.col], cluster = cluster.IDs[feat.inds])
                                        feat.locs
                                    })
        })
    }, chr.blocks.inds, names(chr.blocks.inds)))
}

setMethod("closeClusterFeatures", "ClusteredScoresList", function(c.list, n.consec = NULL, n.perm = 100, name.col = NULL, verbose = TRUE)
{
    if(is.null(n.consec))
        stop("'n.consec' not given.")

    tx.df <- as.data.frame(c.list@anno)
    cluster.IDs <- clusters(c.list)
    pos <- ifelse(tx.df[, "strand"] == '+', tx.df[, "start"], tx.df[, "end"])
    tx.order <- order(tx.df[, "seqnames"], pos)
    tx.df <- tx.df[tx.order, ]
    cluster.IDs <- cluster.IDs[tx.order]
    
    if(verbose)
        message("Finding blocks in experimental clusters.")
    blocks <- .clusterBlocks(tx.df, cluster.IDs, name.col, n.consec)

    if(verbose)
        message("Finding blocks in randomly ordered clusters.")
    rand.blocks <- lapply(1:n.perm, function(x)
                         {
                             rand.IDs <- sample(cluster.IDs, length(cluster.IDs))
                             .clusterBlocks(tx.df, rand.IDs, name.col, n.consec)
                         })
    
    list(blocks = blocks, p.value = length(which(sapply(rand.blocks, length) >= length(blocks))) / n.perm)
})
