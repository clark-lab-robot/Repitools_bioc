setGeneric("colourBlocks", function(blocks, ...){standardGeneric("colourBlocks")})

setMethod("colourBlocks", "GRangesList", function(blocks, cols = NULL, name.col = NULL, verbose = TRUE)
{
    if(is.null(cols))
        stop("Colours to draw features from each cluster in not given.")

    if(length(cols) != length(metadata(blocks)[["clusters"]]))
        stop("Length of colours provided differs to number of clusters used.")

    if(is.null(name.col))
        stop("Name column of 'blocks' meta data not given.")

    up <- metadata(blocks)[["range"]][1]
    down <- metadata(blocks)[["range"]][2]
    max.away <- max(up, down)
    lapply(blocks, function(x)
    {
        x.df <- as.data.frame(x)
        posns <- ifelse(x.df[, "strand"] == '+', x.df[, "start"], x.df[, "end"])
        clusts <- values(x)[, "cluster"]
        x.min <- min(posns) - 2 * max.away
        x.max <- max(posns) + 2 * max.away
        x.mid <- mean(c(x.min, x.max))
        chr  <- x.df[1, "seqnames"]
        plot.new()
        plot.window(c(x.min, x.max), ylim = c(-1, 1))
        axis(1, at = seq(x.min, x.max, length.out = 8))
        text(x.mid, 0, chr)
        rect(rep(x.min, 2), c(-0.25, 0.15), rep(x.max, 2), c(-0.15, 0.25), col = "black")
        text(x.max, -0.3, '-')
        text(x.min, 0.3, '+')
        
        regions.x <- ranges(featureBlocks(x, up, down, keep.strand = TRUE))
        x.lefts <- start(regions.x)
        x.rights <- end(regions.x) 
        regions.y <- as.numeric(ifelse(strand(x) == '+', 0.2, -0.2))
        y.tops <- regions.y + 0.05
        y.bottoms <- regions.y - 0.05
        
        rect(x.lefts, y.bottoms, x.rights, y.tops, col = cols[match(clusts, levels(clusts))])    

        y.lab.pos <- ifelse(x.df[, "strand"] == '+', y.tops + 0.25, y.bottoms - 0.25)
        text(posns, y.lab.pos, x.df[, 5 + name.col], srt = 45)
    })
})
