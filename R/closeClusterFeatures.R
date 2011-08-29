setGeneric("closeClusterFeatures", function(c.list, ...){standardGeneric("closeClusterFeatures")})

.clusterBlocks <- function(tx.df, cluster.IDs, n.consec)
{
    anno.cols <- numeric()
    if(ncol(tx.df) > 5) anno.cols <- 6:ncol(tx.df)
    blocks <- unlist(tapply(cluster.IDs, tx.df[, "seqnames"], function(x)
                    {
                        id.levels <- sort(unique(x))
                        chr.blocks <- lapply(id.levels, function(y)
                        {
                            which.id <- IRanges(as.numeric(names(x[x == y])), width = 1)
                            which.consec <- reduce(which.id)
                            which.consec <- which.consec[width(which.consec) >= n.consec]
                            which.dispersed <- reduce(which.id, min.gapwidth = 2)
                            blocks.inds <- as.list(unique(which.dispersed[subjectHits(findOverlaps(which.consec, which.dispersed))]))
                            if(length(blocks.inds) > 0)
                            {
                                sapply(blocks.inds, function(z)
                                {
                                    cluster.block <- GRanges(tx.df[z, "seqnames"],
                                                             IRanges(tx.df[z, "start"], tx.df[z, "end"]), 
                                                             tx.df[z, "strand"])
                                    values(cluster.block) <- DataFrame(tx.df[z, anno.cols], cluster = cluster.IDs[z])
                                    cluster.block
                                })
                            } else {
                                GRanges()
                            }
                        })
                    }), use.names = FALSE)
    blocks <- blocks[sapply(blocks, function(x) length(x) > 0)]
}

setMethod("closeClusterFeatures", "ClusteredScoresList", function(c.list, n.consec = NULL, n.perm = 100, verbose = TRUE)
{
    if(is.null(n.consec))
        stop("'n.consec' not given.")

    tx.df <- as.data.frame(c.list@anno)
    cluster.IDs <- clusters(c.list)
    pos <- ifelse(tx.df[, "strand"] == '+', tx.df[, "start"], tx.df[, "end"])
    tx.order <- order(tx.df[, "seqnames"], pos)
    tx.df <- tx.df[tx.order, ]
    cluster.IDs <- cluster.IDs[tx.order]
    names(cluster.IDs) <- 1:length(cluster.IDs)
    
    if(verbose)
        message("Finding blocks in experimental clusters.")
    blocks <- GRangesList(.clusterBlocks(tx.df, cluster.IDs, n.consec))
    metadata(blocks) <- list(clusters = levels(cluster.IDs), range = c(c.list@up, c.list@down))

    if(verbose)
        message("Finding blocks in randomly ordered clusters.")
    rand.blocks <- lapply(1:n.perm, function(x)
                         {
                             rand.IDs <- sample(cluster.IDs, length(cluster.IDs))
                             length(.clusterBlocks(tx.df, rand.IDs, n.consec))
                         })
    result <- list(blocks = blocks, p.value = length(which(rand.blocks >= length(blocks))) / n.perm)
    result
})
