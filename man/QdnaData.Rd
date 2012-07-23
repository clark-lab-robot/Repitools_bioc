\name{QdnaData}
\alias{QdnaData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A container for doing ABCD-DNA analyses
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
QdnaData(counts, regions, design, cnv.offsets = NULL, neutral = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{counts}{
%%     ~~Describe \code{counts} here~~
}
  \item{regions}{
%%     ~~Describe \code{regions} here~~
}
  \item{design}{
%%     ~~Describe \code{design} here~~
}
  \item{cnv.offsets}{
%%     ~~Describe \code{cnv.offsets} here~~
}
  \item{neutral}{
%%     ~~Describe \code{neutral} here~~
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
function (counts, regions, design, cnv.offsets = NULL, neutral = NULL) 
{
    require(edgeR)
    stopifnot(!is.null(counts) | !is.null(regions) | !is.null(design))
    stopifnot(length(regions) == nrow(counts))
    stopifnot(nrow(design) == ncol(counts))
    stopifnot(is(regions, "GRanges"))
    if (is.null(cnv.offsets)) 
        cnv.offsets <- matrix(1, nrow = nrow(counts), ncol = ncol(counts))
    if (is.null(neutral)) 
        neutral <- rep(TRUE, nrow(counts))
    stopifnot(all(dim(counts) == dim(cnv.offsets)))
    stopifnot(length(neutral) == nrow(counts))
    suppressMessages(d <- DGEList(counts = counts, genes = as.data.frame(regions)[, 
        -c(4:5)]))
    new("QdnaData", list(DGEList = d, design = as.matrix(design), 
        cnv.offsets = as.matrix(cnv.offsets), neutral = as.logical(neutral), 
        sample.specific.calculated = FALSE))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
