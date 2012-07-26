\name{QdnaData}
\alias{QdnaData}
\alias{QdnaData-class}
\alias{class:QdnaData}
\alias{show,QdnaData-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A container for quantitative DNA sequencing data for ABCD-DNA analyses
}
\description{
\code{QdnaData} objects form the basis for differential analyses of quantitative DNA sequencing data(i.e. ABCD-DNA).  A user is required to specify the minimum elements: a count table, a list of regions and a design matrix.  For copy-number-aware analyses, a table of offsets and the set of neutral regions needs to be given.
}
\usage{
QdnaData(counts, regions, design, cnv.offsets = NULL, neutral = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{counts}{
table of counts for regions of interest across all samples
}
  \item{regions}{
a \code{GRanges} object giving the regions
}
  \item{design}{
a design matrix
}
  \item{cnv.offsets}{
a table of offsets.  If unspecified (or \code{NULL}), a matrix of 1s (i.e. no CNV) is used
}
  \item{neutral}{
a logical vector, or indices, of the regions deemed to be neutral. If unspecified (or \code{NULL}), all regions are used
}
}
\details{
\code{QdnaData} objects are geared for general differential analyses of qDNA-seq data.  If CNV is present and prominent, the objects and methods available with \code{QdnaData} perform adjustments and spot checks before the differential analysis.
}
\value{
a \code{\link{QdnaData}} object (effectively a list) is returned
}
\references{
http://imlspenticton.uzh.ch/robinson_lab/ABCD-DNA/ABCD-DNA.html
}
\author{
Mark Robinson
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{getSampleOffsets}}, \code{\link{plotQdnaByCN}}, \code{\link{setCNVOffsets}}
}
\examples{
require(GenomicRanges)
cnt <- matrix(rpois(20,lambda=10),ncol=4)
gr <- GRanges("chr1",IRanges(seq(2e3,6e3,by=1e3), width=500))
des <- model.matrix(~c(0,0,1,1))
qd <- QdnaData( counts=cnt, regions=gr, design=des)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
