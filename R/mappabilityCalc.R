setGeneric("mappabilityCalc", function(x, organism, ...){standardGeneric("mappabilityCalc")})

.mappabilityBSGenome <- function(regions.by.chr, organism)
{
    chr.maxs <- seqlengths(organism)[names(regions.by.chr)]
    
    mappability.by.chr <- mapply(function(y, z)
    {
        # Handle case of windows overlapping past ends of chromosome.
        inside.regions <- restrict(y, 1, z, keep.all.ranges = TRUE)
        window.seqs <- suppressWarnings(getSeq(organism, inside.regions))
        unmap.counts <- letterFrequency(DNAStringSet(window.seqs), 'N')
        unmap.counts <- unmap.counts + width(y) - width(inside.regions)
        1 - (unmap.counts / width(y))
    }, as.list(regions.by.chr), as.list(chr.maxs), SIMPLIFY = FALSE)
}

.mappabilityFASTA <- function(regions.by.chr, map.file, verbose)
{
    file.conn <- file(map.file, 'r')
    chr.found <- FALSE
    chr.map <- list()

    repeat
    {
        text.line <- readLines(file.conn, 1)
        if(grepl("^~[^~]", text.line) && chr.found == FALSE)
        # First time seeing a chromosome.
        {
            chr.name <- gsub('~', '', text.line)
            chr.found <- TRUE
            map.string <- character(10000000)
            fileLineIndex = 1
            vectorIndex = 1
            if(verbose) message("Reading mappability for ", chr.name, '.')
        } else if((grepl("^~[^~]", text.line) && chr.found == TRUE) || length(text.line) == 0) 
        {
        # End of reading data for a chromosome. Process user regions.
            if(verbose) message("Finished reading mappability for ", chr.name, '.')
            map.string <- paste(map.string, collapse = '')
	    if(chr.name %in% names(regions.by.chr))
            {
                chr.regions <- regions.by.chr[[chr.name]]
                inside.regions <- restrict(chr.regions, 1, nchar(map.string),
                                           keep.all.ranges = TRUE)
                map.set <- BStringSet(map.string, start(inside.regions),
                                      width(inside.regions))
                unique.bases <- letterFrequency(map.set, '!')
                chr.map <- c(chr.map, list(unique.bases / width(chr.regions)))
                names(chr.map)[length(chr.map)] <- chr.name
                if(verbose) message("Processed regions on ", chr.name, '.')
            }
            if(length(text.line) == 0) break 
            map.string <- character(10000000)
            vectorIndex <- 1
            chr.name <- gsub('~', '', text.line)
            if(verbose) message("Reading mappability for ", chr.name, '.')
        } else if (chr.found == TRUE) {
        # Line of mappability scores. Append to existing scores list.
            map.string[vectorIndex] <- text.line
            vectorIndex <- vectorIndex + 1
            fileLineIndex <- fileLineIndex + 1
        }
    }

    close(file.conn)
    if(length(setdiff(names(regions.by.chr), names(chr.map))) > 0)
        stop("Some user regions' chromosomes are not present in the FASTA file.")
    chr.map[names(regions.by.chr)]
}
    
setMethod("mappabilityCalc", c("GRanges", "MappabilitySource"),
          function(x, organism, window = NULL, type = c("block", "TSS", "center"),
          verbose = TRUE)
{
    type <- match.arg(type)
    if(type == "block" && !is.null(window))
        stop("'window' is meaningless when region type is \"block\".")

    if(type %in% c("TSS", "center") && is.null(window))
        stop("Using a reference point but window size around it was not specified.")

    info <- "Calculating mappability"
    if(type == "block") info <- paste(info, "in supplied blocks.")
    if(type == "TSS") info <- paste(info, window, "bases around TSSs.")
    if(type == "center") info <- paste(info, window, "bases around feature centres.")
    
    if(verbose) message(info)

    if(type %in% c("TSS", "center"))
    {
        if(type == "TSS")
            x.posns <- IRanges(as.numeric(ifelse(strand(x) == '+', start(x), end(x))),
                           width = 1)
        if(type == "center")
            x.posns <- IRanges(as.integer((start(x) + end(x)) / 2), width = 1)
        
        ranges(x) <- x.posns
        x <- resize(x, window, "center")
    }

    strand(x) <- "+"
    chrs <- as.character(seqnames(x))
    regions.by.chr <- split(x, chrs)

    if(class(organism) == "character")
    {
        unsplit(.mappabilityFASTA(regions.by.chr, organism, verbose), chrs)
    } else {
        unsplit(.mappabilityBSGenome(regions.by.chr, organism), chrs)
    }
})
    
setMethod("mappabilityCalc", c("data.frame", "MappabilitySource"),
         function(x, organism, window = NULL, type = c("block", "TSS", "center"), ...)
{
    type <- match.arg(type)
    if(type == "block" && !is.null(window))
        stop("'window' is meaningless when region type is \"block\".")
    
    if(type %in% c("TSS", "center") && is.null(window))
        stop("Using a reference point but window size around it was not specified.")

    if(type == "block")
    {
        x <- GRanges(x$chr, IRanges(x$start, x$end))
    } else {
        if (is.null(x$position))
        {
            if(type == "TSS")
                x$position <- ifelse(x$strand == '+', x$start, x$end)
            if(type == "center")
                x$position <- as.integer((x$start + x$end) / 2)
        }
        x <- GRanges(x$chr, IRanges(x$position, width = 1))
        x <- resize(x, window, "center")
    }

    mappabilityCalc(x, organism, window = NULL, type = "block", ...)
})
