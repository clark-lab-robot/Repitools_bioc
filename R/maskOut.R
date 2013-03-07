maskOut <- function(x, ranges){

    if(class(ranges) != "GRanges"){
        stop("ranges must be of class GRanges")
    }

    fo <- findOverlaps(ranges, windows(x))

    mask <- logical(length(x))
    mask[fo@subjectHits] <- TRUE

    maskEmpBayes(x) <- mask
    return(x)
}
