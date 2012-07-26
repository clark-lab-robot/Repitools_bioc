\name{getSampleOffsets}
\alias{getSampleOffsets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the sample-specific offsets, using the neutral state
}
\description{
ABCD-DNA combines CNV offsets with sample specific factors.  This function calculates the latter, using a set of neutral regions (and corresponding counts in the count table).
}
\usage{
getSampleOffsets(obj, ref = 1, quantile = 0.99, min.n = 100, plot.it = FALSE, force = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
a \code{\link{QdnaData}} object
}
  \item{ref}{
integer index, giving the sample to use as reference
}
  \item{quantile}{
  quantile of the A-values to use
}
  \item{min.n}{
  minimum number of points to include
}
  \item{plot.it}{
  logical, whether to plot an M-A plot for each sample against the reference (default: \code{FALSE})
}
  \item{force}{
  logical, whether to recalculate the sample-specific offsets (only needed if they are already calculated)
}
  \item{\dots}{
  arguments to pass to the \code{maPlot} function
}
}
\details{
  The sample-specific offset is calculated as the median M-value beyond (i.e. to the right) an A-value quantile, using only the copy-number-neutral regions, as specified in the incoming \code{QdnaData} object.
}
\value{
  returns a \code{QdnaData} object (copied from the \code{obj} argument) and populates the \code{$DGEList$samples$norm.factors} element and sets the \code{$sample.specific.calculated} to \code{TRUE}.
}
\references{
http://imlspenticton.uzh.ch/robinson_lab/ABCD-DNA/ABCD-DNA.html
}
\author{
Mark Robinson
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{QdnaData}}
}
\examples{
# library(Repitools)
# qd <- QdnaData(counts=counts, regions=gb, design=design, 
#               cnv.offsets=cn, neutral=(regs=="L=4 P=2"))
# qd <- getSampleOffsets(qd,ref=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
