setGeneric("absoluteCN", function(input.windows, ip.windows, ...){standardGeneric("absoluteCN")})

setMethod("absoluteCN", c("GRanges", "GRanges"),
    function(input.windows, ip.windows, input.counts = NULL, gc.params = NULL, verbose = TRUE)
{
    if(is.null(gc.params))
        stop("gc.params must be specified for absolute copy estimation.")

    n.bins <- gc.params@n.bins

    # Find which counting windows have sufficient mappability.
    input.win.mappability <- mappabilityCalc(input.windows, gc.params@mappability) * 100
    mappable <- input.win.mappability > gc.params@min.mappability
    input.windows <- input.windows[mappable, ]
    input.counts <- input.counts[mappable, ]
    input.counts <- input.counts * 100 / input.win.mappability[mappable]
    abs.CN <- apply(input.counts, 2, function(x) x / median(x))

    # Get the GC content of windows.
    gc <- gcContentCalc(input.windows, gc.params@genome)

    # Break GC content into bins, and find mode of bins, also including adjacent bins.
    gc.range <- range(gc)
    bins <- seq(gc.range[1], gc.range[2], length.out = n.bins) # Midpoint of each bin.
    mode.gcs <- apply(abs.CN, 2, function(x)
    {
        modes <- rep(NA, n.bins)
        for(index in 2:(gc.params@n.bins-1))
        {
            in.bins <- which(gc >= bins[index - 1] & gc < bins[index + 1])
            count.dens <- density(x[in.bins])
            modes[index] <- as.numeric(count.dens$x[which.max(count.dens$y)])
        }
        modes
    })

    # Find expected copy number for each window, based on its GC content.
    model.CN <- apply(mode.gcs, 2, function(x)
    {
        counts.model <- lm(x ~ poly(bins, gc.params@poly.degree))
        predict(counts.model, data.frame(bins = gc))
    })

    # Adjust the real counts by dividing by expected counts.
    scaled.CN <- t(t(abs.CN / model.CN) * gc.params@ploidy)
    
    if(verbose == TRUE) message("Mapping copy number windows to IP windows.")
    map <- findOverlaps(ip.windows, input.windows, select = "first")

    features.CN <- scaled.CN[map, ]
    rownames(features.CN) <- .getNames(ip.windows)

    features.CN
})

setMethod("absoluteCN", c("data.frame", "data.frame"),
    function(input.windows, ip.windows, ...)
{
    absoluteCN(annoDF2GR(input.windows), annoDF2GR(ip.windows), ...)
})
