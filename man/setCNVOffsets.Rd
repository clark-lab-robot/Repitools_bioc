\name{setCNVOffsets}
\alias{setCNVOffsets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set the CNVOffsets of a \code{QdnaData} object
}
\description{
A utility function to manually add CNV offset to a \code{QdnaData} object
}
\usage{
setCNVOffsets(obj, cnv.offsets)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
a \code{QdnaData} object
}
  \item{cnv.offsets}{
a \code{matrix} of offsets (presumably copy number)
}
}
\value{
a \code{QdnaData} object
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
#                neutral=(regs=="L=4 P=2"))
# qd <- setCNVoffsets(qd, cn)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
