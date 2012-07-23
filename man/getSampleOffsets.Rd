\name{getSampleOffsets}
\alias{getSampleOffsets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the sample-specific offsets, using the neutral state
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getSampleOffsets(obj, ref = 1, quantile = 0.99, min.n = 100, plot.it = FALSE, force = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
%%     ~~Describe \code{obj} here~~
}
  \item{ref}{
%%     ~~Describe \code{ref} here~~
}
  \item{quantile}{
%%     ~~Describe \code{quantile} here~~
}
  \item{min.n}{
%%     ~~Describe \code{min.n} here~~
}
  \item{plot.it}{
%%     ~~Describe \code{plot.it} here~~
}
  \item{force}{
%%     ~~Describe \code{force} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (obj, ref = 1, quantile = 0.99, min.n = 100, plot.it = FALSE, 
    force = FALSE, ...) 
{
    if (!is(obj, "QdnaData")) 
        stop("Only operates on 'QdnaData' objects.")
    if (obj$sample.specific.calculated & !force) 
        stop("Sample-specific factors have already been calculated. Set force=TRUE to recalculate.")
    cnt <- obj$DGEList$counts
    nc <- ncol(cnt)
    nf <- rep(0, nc)
    idx <- setdiff(1:nc, ref)
    rs <- obj$DGEList$samples$lib.size
    neu <- obj$neutral
    for (i in idx) {
        map <- maPlot(cnt[neu, ref], cnt[neu, i], normalize = TRUE, 
            plot.it = plot.it, ...)
        q <- .getQuantile(map, min.n, quantile)
        if (plot.it) {
            grid()
            abline(v = q, col = "blue")
        }
        nf[i] <- median(map$M[map$A > q])
        if (plot.it) 
            abline(h = nf[i], col = "red", lwd = 4)
    }
    obj$DGEList$samples$norm.factors <- exp(nf - mean(nf))
    obj$sample.specific.calculated <- TRUE
    obj
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
