setGeneric("summarizeScores", function(scores.list, design, ...)
                                 {standardGeneric("summarizeScores")})

setMethod("summarizeScores", c("ScoresList", "matrix"),
    function(scores.list, design, verbose = TRUE)
{

    list.scores <- tables(scores.list)
    n.samples <- length(list.scores)
    if(n.samples != nrow(design))
        stop("Design matrix has different number of rows to number of samples",
             " in scores list.")

    if(!is.null(scores.list@s.width))
    {
        smoothings <- apply(design, 2, function(x)
                           {
                               smooth <- unique(scores.list@s.width[x != 0])
                               if(length(smooth) > 1)
                                   stop("Trying to subtract scores with samples ",
                                        "that have different smoothing widths.")
                               smooth
                           })

    } else {
        smoothings <- NULL
    }

    weights.matrix <- apply(design, 2, function(x)
                           {
                               x[x == 1] <- 1 / sum(x == 1)
                               x[x == -1] <- -1 / sum(x == -1)
                               x
                           })



    if(verbose) message("Calculating feature scores differences.")
    
    n.rows <- nrow(list.scores[[1]])
    n.cols <- ncol(list.scores[[1]])
    wts.list <- split(weights.matrix, rep(1:ncol(weights.matrix), each = nrow(weights.matrix)))
    summaries <- lapply(wts.list, function(x)
                 {
                     s.matrix <- matrix(0, nrow = n.rows, ncol = n.cols)
                     for(index in 1:n.samples)
                     {
                        s.matrix <- s.matrix + list.scores[[index]] * x[index] 
                     }
                     s.matrix
                 })

    scores.list@names <- colnames(design)
    scores.list@scores <- summaries
    scores.list@s.width <- smoothings

    scores.list
})
