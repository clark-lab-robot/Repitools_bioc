setGeneric("GCadjustCopy", function(input.windows, input.counts, gc.params, ...)
                                 {standardGeneric("GCadjustCopy")})

setMethod("GCadjustCopy", c("GRanges", "matrix", "GCAdjustParams"),
    function(input.windows, input.counts, gc.params, verbose = TRUE)
{
    n.bins <- gc.params@n.bins

    if(length(input.windows) != nrow(input.counts))
        stop("Number of input windows different to rows of counts.")

    # Find which counting windows have sufficient mappability.
    input.win.mappability <- mappabilityCalc(input.windows, gc.params@mappability)
    mappable <- input.win.mappability * 100 > gc.params@min.mappability 
    counts.map.scaled <- input.counts / input.win.mappability
    abs.CN <- apply(counts.map.scaled, 2, function(x) x / median(x[mappable]))

    # Get the GC content of windows.
    gc <- gcContentCalc(input.windows, gc.params@genome)

    # Break GC content into bins, and find mode of bins, also including adjacent bins.
    gc.range <- range(gc[mappable])
    bins <- seq(gc.range[1], gc.range[2], length.out = n.bins) # Midpoint of each bin.
    mode.gcs <- apply(abs.CN, 2, function(x)
    {
        modes <- rep(NA, n.bins)
        for(index in 2:(gc.params@n.bins-1))
        {
            in.bins <- which(gc >= bins[index - 1] & gc < bins[index + 1] & mappable)
            if(length(in.bins) >= gc.params@min.bin.size)
            {
                count.dens <- density(x[in.bins])
                modes[index] <- as.numeric(count.dens$x[which.max(count.dens$y)])
            }
        }
        modes
    })

    # Find polynomial models for each sample.
    if(verbose) message("Fitting models to counts in GC bins.")
    models <- apply(mode.gcs, 2, function(x)
    {
        lm(x ~ poly(bins, gc.params@poly.degree))
    })
    names(models) <- colnames(input.counts)

    # Estimate single CN score for each sample for each window.
    if(verbose) message("Estimating single copy score with each model for all windows.")
    single.CNs <- sapply(models, function(x) predict(x, data.frame(bins = gc)))

    # Adjust the real counts by dividing by expected counts.
    adj.CN <- t(t(abs.CN / single.CNs) * gc.params@ploidy)
    adj.CN[!mappable, ] <- NA

    CopyEstimate(gc.params@ploidy, input.windows, input.counts, input.win.mappability,
                 gc, models, adj.CN)
})

setMethod("GCadjustCopy", c("data.frame", "matrix", "GCAdjustParams"),
    function(input.windows, input.counts, gc.params, ...)
{
    GCadjustCopy(annoDF2GR(input.windows), input.counts, gc.params, ...)
})
