\name{maskOut}
\alias{maskOut}
\title{
    Function to mask suspicious regions.
}
\description{
    Function to mask out regions that are prone to causing problems
    in the empirical Bayes approach \code{empBayes}. 
    The corresponding bins are marked and in the empirical Bayes 
    approach not taken into account. Notice that methylation estimates 
    using \code{methylEst} will nevertheless be produced 
    for these bins.
}
\usage{
maskOut(x, ranges)
}
\arguments{
    \item{x}{
    Object of class \code{BayMethList}.
}
    \item{ranges}{
    A GRanges object definining the coordinates of regions to be masked out.
}
}
\value{
A \code{BayMethList} object where the slot \code{maskout} is filled with a
boolean vector indicating which bins will be excluded in \code{empBayes}.
}
\author{
    Andrea Riebler
}
\examples{
    if(require(BSgenome.Hsapiens.UCSC.hg18)){
        windows <- genomeBlocks(Hsapiens, chrs="chr21", width=100, spacing=100)
        cpgdens <- cpgDensityCalc(windows, organism=Hsapiens, 
            w.function="linear", window=700)  
        co <- matrix(rnbinom(length(windows), mu=10, size=2), ncol=1)
        sI <- matrix(rnbinom(2*length(windows), mu=5, size=2), ncol=2)
        bm <- BayMethList(windows=windows, control=co, 
            sampleInterest=sI, cpgDens=cpgdens)

        # mask out unannotated high copy number regions
        # see Pickrell et al. (2011), Bioinformatics 27: 2144-2146.
        data(hcRegions)

        bm <- maskOut(bm, hcRegions)
   }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{programming}
