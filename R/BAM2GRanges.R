setGeneric("BAM2GRanges", signature = "path", function(path, ...)
{standardGeneric("BAM2GRanges")})
setGeneric("BAM2GRangesList", signature = "paths", function(paths, ...)
{standardGeneric("BAM2GRangesList")})

setMethod("BAM2GRanges", "character",
          function(path, what = character(),
                   flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE),
                   verbose = TRUE)
          {
              if(length(path) > 1)
                  stop("This method is only for one BAM file. See ?BAM2GRangesList.")
              if(!grepl("bam", path))
                  stop("Path does not seem to be for a BAM file.")
            
              if(verbose == TRUE)
                  cat("Reading BAM file ", basename(path), ".\n", sep = '')
              filters <- ScanBamParam(what = what, flag = flag)
	      aligns <- readGAlignments(path, param = filters)
	      aligns.info <- values(aligns)
              aligns.gr <- as(aligns, "GRanges")
	      values(aligns.gr) <- aligns.info
	      aligns.gr
          })

setMethod("BAM2GRangesList", "character",
          function(paths, what = character(),
                   flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE),
                   verbose = TRUE)
          {
              GRangesList(lapply(paths, function(x) BAM2GRanges(x, what, flag, verbose)))
          })
