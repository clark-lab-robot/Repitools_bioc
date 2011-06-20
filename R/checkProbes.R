setGeneric("checkProbes", signature = c("regs", "probes"), function(regs, probes, ...)
                                           {standardGeneric("checkProbes")})

setMethod("checkProbes", c("GRanges", "GRanges"),
	function(regs, probes, up = NULL, down = NULL, ...)
{
    require(GenomicRanges)

    str.regs <- as.character(strand(regs))
    chr.regs <- as.character(seqnames(regs))
    pr.names <- elementMetadata(probes)[, "name"]
    reg.names <- elementMetadata(regs)[, "name"]

    ref.pos <- ifelse(str.regs == '+', start(regs), end(regs))
    end.pos <- ifelse(str.regs == '+', end(regs), start(regs))

    if(!is.null(up) & !is.null(down))
    {
	start(regs) <- ifelse(str.regs == '+', ref.pos - up, ref.pos - down)
	start(regs)[start(regs) < 0] <- 0
	end(regs) <- ifelse(str.regs == '+', ref.pos + down, ref.pos + up)
    }

    hits <- table(pr.names)
    strand(regs) <- '*'

    map <- findOverlaps(regs, probes)@matchMatrix
    inds <- split(map[, 2], factor(map[, 1], levels = 1:length(regs)))

    invisible(mapply(function(x, y) {
	par(oma = c(4, 1, 1, 1))
	pr.hits <- match(pr.names[y], names(hits))
        if(length(pr.hits) > 0)
        {
	    x.pos <- split(c(start(probes[y]), end(probes[y])),
                           rep(1:(length(probes[y])), 2))
	    y.max = max(pr.hits)

            # Pretty printing of y-axis.
            if(y.max > 100)
	    {
		ticks = c(1, seq(25, y.max, 25))
	    } else if(y.max > 10)
	    {
		ticks = c(1, seq(5, y.max, 5))
	    } else {
		ticks = 1:y.max
	    }
            x.lim <- c(start(regs)[x], end(regs)[x])
	    matplot(x.pos[[1]], rep(pr.hits[1], 2), ylim = c(0, y.max),
                    type = 'l', xlim = x.lim, yaxt = 'n',
                    main = reg.names[x], xlab = "Position",
                    ylab = "Probe Hits", ...)
	    axis(2, ticks, las = 2)
            if(length(y) > 1)
            {
		mapply(function(v, w)
		{
		    matlines(v, rep(w, 2), ...)
		}, x.pos[2:length(x.pos)], pr.hits[2:length(pr.hits)])
            }
        } else { # Plot empty plot for unprobed gene.
                r.st <- start(regs[x])
                r.end <- end(regs[x])
                plot(1, type = 'n', main = reg.names[x], xlim = c(r.st, r.end),
                     ylim = c(0, 1), yaxt = 'n', xlab = "Position", ylab = "Probe Hits")
                axis(2, 0:1, las = 2)
                text(ifelse(is.null(up), (r.st + r.end) / 2, 1.1 * r.st),
                         0.5, "No Hits", col = "red")
        }
	
	mtext("TSS", 3, 0, at = ref.pos[x], col = "red")
	abline(v = ref.pos[x], col = "red")
	mtext("TTS", 3, 0, at = end.pos[x], col = "red")
	abline(v = end.pos[x], col = "red")
	abline(h = 1, lty = 2, col = "darkgrey")
	mtext(chr.regs[x], 1, font = 2, at = 0.52, outer = TRUE)
    }, as.numeric(names(inds)), inds))
})

setMethod("checkProbes", c("data.frame", "data.frame"),
    function(regs, probes, up = NULL, down = NULL, ...)
{
    checkProbes(annoDF2GR(regs), annoDF2GR(probes), up = up, down = down, ...)
})
