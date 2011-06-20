setGeneric("featureBlocks", signature = "anno", function(anno, ...)
                                {standardGeneric("featureBlocks")})

setMethod("featureBlocks", "GRanges",
    function(anno, up = NULL, down = NULL, dist = c("base", "percent"), keep.strand = FALSE)
{
    require(GenomicRanges)
    dist <- match.arg(dist)
    .validate(anno, up, down)

    str <- as.character(strand(anno))
    f.st <- start(anno)
    f.end <- end(anno)
    wd <- width(anno)
    pos <- str == '+'

    if(all(str == '*'))
        ref.points <- round((f.st + f.end) / 2)
    else
        ref.points <- ifelse(pos, f.st, f.end)

    if(dist == "percent")
    {
        starts = as.numeric(ifelse(pos, ref.points - wd * up/100,
                                        ref.points - wd * down/100))
        ends = as.numeric(ifelse(pos, ref.points + wd * down/100,
                                        ref.points + wd * up/100))
    } else {
	starts = as.numeric(ifelse(pos, ref.points - up,
                                        ref.points - down))
	ends = as.numeric(ifelse(pos, ref.points + down,
                                      ref.points + up))
    }
    f.names <- .getNames(anno)

    GRanges(seqnames(anno), IRanges(starts, ends), if(keep.strand) str else '*',
            name = f.names)
})

setMethod("featureBlocks", "data.frame",
    function(anno, ...)
{
    featureBlocks(annoDF2GR(anno), ...)
})
