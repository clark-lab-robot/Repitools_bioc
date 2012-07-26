\name{abcdDNA}
\alias{abcdDNA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A wrapper for fitting the offset-adjusted ABCD-DNA GLM
}
\description{
This function performs differential analyses, given a \code{QdnaData} object with the sample-specific offsets already calculated (i.e. call \code{getSampleOffsets} before calling \code{abcdDNA}), a coefficient (or set of coefficients) to test and dispersion(s).  In essence, the function is a wrapper for constructing the offset matrix, fitting the generalized linear model and performing a likelihood ratio test.
}
\usage{
abcdDNA(obj, coef = ncol(obj$design), dispersion = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
a \code{QdnaData} object
}
  \item{coef}{
coefficient (or coefficients) of the design matrix to test
}
  \item{dispersion}{
estimate(s) of dispersion to use for negative binomial testing
}
}
\details{
This function is simply a wrapper for taking the details in an \code{QdnaData} object and perform the differential analyses, adjusting for copy number if specified.
}
\value{
a \code{DGEGLM} (see the edgeR package) containing the results of the differential comparison
}
\references{
http://imlspenticton.uzh.ch/robinson_lab/ABCD-DNA/ABCD-DNA.html
}
\author{
Mark Robinson
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{QdnaData}}, 
}
\examples{
# library(Repitools)
# qd <- QdnaData(counts=counts, regions=gb, design=design, 
#                cnv.offsets=cn, neutral=(regs=="L=4 P=2"))
# qd <- getSampleOffsets(qd,ref=1)
# plotQdnaByCN(qd,cnv.group=regs,idx.ref=3,idx.sam=2)
# f <- abcdDNA(qd, dispersion=.05, coef=2)
# topTags(f)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
