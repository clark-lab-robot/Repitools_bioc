setGeneric("genomeBlocks", function(genome, ...) {standardGeneric("genomeBlocks")})

setMethod("genomeBlocks", "numeric",
    function(genome, chrs = names(genome), width = NULL, spacing = width)
{
    if(is.null(width))
        stop("Block width not given.")

    require(GenomicRanges)
    
    chr.windows <- lapply(chrs, function(x) {
      centres <- seq.int(min(spacing/2,genome[x]),genome[x],spacing)
      GRanges(seqnames = names(genome[x]),
              ranges = IRanges(start = pmax(centres-width/2+1,1), 
                               end=pmin(centres+width/2, genome[x])))
    })
    
    suppressWarnings(do.call(c, chr.windows))
})

setMethod("genomeBlocks", "BSgenome",
    function(genome, chrs = seqnames(genome), width = NULL, spacing = width)
{
    if(is.null(width))
        stop("Block width not given.")

    require(BSgenome)
    
    chr.lengths <- seqlengths(genome)[chrs]
    genomeBlocks(chr.lengths, chrs = names(chr.lengths), width = width, spacing = spacing)
})
