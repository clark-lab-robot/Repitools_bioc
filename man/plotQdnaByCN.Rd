\name{plotQdnaByCN}
\alias{plotQdnaByCN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting the response of qDNA-seq data by CNV
}
\description{
Given groupings of relative CNV state, this function produces M-A (log-fold-change versus log-average) plots to compare two samples relative read densities.  In addition, it calculates a scaling factor at a specified quantile and plots the median M value across all the groups.
}
\usage{
plotQdnaByCN(obj, cnv.group, idx.ref = 1, idx.sam = 2, min.n = 100, quantile = 0.99, ylim = c(-5, 5), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
a \code{QdnaData} object
}
  \item{cnv.group}{
a \code{character} vector or \code{factor} giving the relative CNV state.  This must be the same length as the number of regions in \code{obj}
}
  \item{idx.ref}{
index of the reference sample (denominator in the calculation of M values)
}
  \item{idx.sam}{
index of the sample of interest (numerator in the calculation of M values)
}
  \item{min.n}{
  minimum number of points to include
}
  \item{quantile}{
  quantile of the A-values to use
}
  \item{ylim}{
y-axis limits to impose on all M-A plots
}
  \item{\dots}{
further arguments sent to \code{maPlot}
}
}
\value{
a plot to the current graphics device
}
\references{
http://imlspenticton.uzh.ch/robinson_lab/ABCD-DNA/ABCD-DNA.html
}
\author{
Mark Robinson
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{QdnaData}}, ~~~
}
\examples{
# library(Repitools)
# qd <- QdnaData(counts=counts, regions=gb, design=design, 
#                cnv.offsets=cn, neutral=(regs=="L=4 P=2"))
# plotQdnaByCN(qd,cnv.group=regs,idx.ref=3,idx.sam=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
