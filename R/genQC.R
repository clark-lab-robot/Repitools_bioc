setGeneric("genQC", function(qc.data, ...) {standardGeneric("genQC")})

.plotFreqs <- function(freqs, l.names, aligned)
{
  layout(matrix(c(1:8), ncol = 4, byrow = TRUE))
    invisible(mapply(function(y, z){
        matplot(y = y, main = paste(ifelse(aligned == TRUE, "Aligned", "All"),
                                    "Bases For", z),
                xlab = "Cycle", ylab = "Base", type = "l", lty = 1, lwd = 2,
                xlim = c(1, nrow(y)), ylim = c(0.00, 0.40))
        abline(h = 0.25, col = "red")
        legend("topright", legend = c("G", "A", "T", "C", "N"), lty = 1, lwd = 2,
               col = 1:5, cex = 0.5)
    }, freqs, l.names))
}

setMethod("genQC", "SequenceQCSet", 
    function(qc.data, expt = "Experiment")
{
    if(is.null(names(qc.data)))
        l.names <- paste("Lane", 1:length(qc.data))
    else
        l.names <- names(qc.data)

    n.reads <- lapply(qc.data, function(x)
                   as.integer(Basic_Statistics(x, "Unaligned")[3, 2]))
    n.aln <- lapply(qc.data, function(x)
                 as.integer(Basic_Statistics(x, "Aligned")[3, 2]))

    avg.q <- sapply(qc.data, function(x)
                 Per_base_sequence_quality(x, "Unaligned")[, 2])

    base.freq <- mapply(function(x, y)
    {
        cbind(Per_base_sequence_content(x, "Unaligned")[, 2:5] / y,
              Per_base_N_content(x, "Unaligned")[, 2] / 100)
    }, qc.data, n.reads, SIMPLIFY = FALSE)

    base.freq.aln <- mapply(function(x, y)
    {
        cbind(Per_base_sequence_content(x, "Aligned")[, 2:5] / y,
              Per_base_N_content(x, "Aligned")[, 2] / 100)
    }, qc.data, n.aln, SIMPLIFY = FALSE)

    mm.freqs <- lapply(qc.data, function(x)
    {
	mism <- Mismatches(x)[, 2:17]
        t(apply(mism, 1, function(y) y / sum(y)))
    })

    dupls <- lapply(qc.data, function(x)
    {
	dupl <- Sequence_Duplication_Levels(x)
	return(as.numeric(dupl[2:nrow(dupl), 2]))
    })
    
    offset <- max(avg.q * 0.1)
    y.lim <- c(min(avg.q) - offset, max(avg.q + offset))
    matplot(avg.q, type = "l", lty = 1, lwd = 2, col = 1:8,
        main = paste("Average Quality Over Cycles for", expt),
        xlab = "Cycle", ylab = "Quality", ylim = y.lim)
        legend("topright", legend = l.names, lty = 1, lwd = 2,
               col = 1:length(l.names))
    .plotFreqs(base.freq, l.names, FALSE)
    .plotFreqs(base.freq.aln, l.names, TRUE)
    par(oma = c(1, 1, 1, 1))
    layout(matrix(c(1:4), ncol = 2, byrow = TRUE), widths=c(5,1))
    invisible(mapply(function(y, z){
        par(mai = c(1,1.2,0.2,0.5))
        cols.mm <- c("black", "darkblue", "purple", "darkgreen", "cyan",
                     "lightgreen", "red", "yellow")
        matplot(y = y,  main = paste("Mismatches For", z), xlab = "Cycle",
                ylab = "Percent of All Mismatches At Cycle", type = "l",
                lty = rep(c(1, 2), each = 8), lwd = 2, xlim = c(1, nrow(y)),
                ylim = c(0.00, 1.00), col = cols.mm)
        abline(v = 1:nrow(y), col = c("lightgrey", "darkgrey"))
        par(mai = c(0,0,0.2,0))
        plot.new()
        legend("topleft", title = "Reference > Read", legend = colnames(y),
               lty = rep(c(1, 2), each = 8), lwd = 2, col = cols.mm)
    }, mm.freqs, l.names))

    layout(matrix(c(1:4), ncol = 2, byrow = TRUE))
    invisible(mapply(function(y, z){
        par(mai = c(0.5,0.5,0.5,0.5))
        bar.x <- barplot(y,  main = paste("Duplication For", z),
                        xlab = "Occurrences", ylab = "Relative Number", xaxt = 'n')
        axis(1, at = bar.x, labels = c(1:9, "10+"), cex.axis = 0.5)
    }, dupls, l.names))

})

setMethod("genQC", "character", 
    function(qc.data, ...)
{
    qual <- lapply(qc.data, function(x) get(load(x)))

    QCset <- new("SequenceQCSet", qual)
    names(QCset) <- names(qc.data)
    genQC(QCset, ...)
})
