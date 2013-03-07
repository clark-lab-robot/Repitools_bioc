\name{methylEst}
\alias{methylEst}
\title{
  Function to derive regional methylation estimates.
}
\description{
Posterior mean and variance for the regional methylation level are derived 
for all genomic regions. Credible intervals can be computed either numerically 
from the posterior marginal distribution or by computing them on logit
scale and transferring them back.
}
\usage{
methylEst(x, controlCI = list(compute = FALSE, method = "Wald", 
    level = 0.95, nmarg = 512, ncpu = NULL))
}
\arguments{
  \item{x}{
    Object of class \code{BayMethList}.
    }
  \item{controlCI}{
    list defining whether credible intervals should be derived. 
    \describe{
    \item{-}{
        \code{compute}
        logical. If 'TRUE' credible intervals are derived. (default FALSE)
    }
    \item{-}{
        \code{method}
        'Wald' (default) for Wald-type intervalse computed on logit scale 
        and transferred back. Two more credible types will come: 'HPD' for 
        highest posterior density intervals; 'quantile' for quantile-based
        intervals;
    }
    \item{-}{
        \code{level}
        numerical value defining the credible level. Default: 0.95.
        }
    \item{-}{
        \code{nmarg}
        Number of points at which the posterior marginal is evaluated
        (only relevant for method="quantile" or method="HPD")
        TODO: Ooops not good!
    }
    \item{-}{
        \code{ncpu}
        Number of CPUs on your machine you would like to use in parallel. 
        If \code{ncpu} is set to NULL, half of the CPUs will be used on
        machines with a maximum of four CPUs, and 2/3 will be used if
        more are available.
        }
    }  
}
}
\details{
The posterior mean and the variance are analytically available and therefore 
straightforward to efficiently compute; credible intervals can be computed
numerically from the posterior marginal distribution.
}
\value{
    A \code{BayMethList} object where the slot \code{methEst} is filled with a
    list containing the following elements:
    \item{mean }{Matrix where the number of columns equals the number of
    samples of interest. Each column contains the posterior mean methylation
    level for each bin.}
    \item{var }{Matrix where the number of columns equals the number of samples 
        of interest. Each column contains posterior variance for each bin.}
    \item{ci }{Credible quantiles}
    \item{W }{Normalisation factor of the posterior marginal distribution 
        for each bin.}
    \item{al }{Prior shape parameter for each bin}
    \item{bl }{Prior scale parameter for each bin}
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

        bm <- determineOffset(bm)
        # should take about 3 minutes for both samples of interests with 2 CPUs.
        bm <- empBayes(bm) 
        bm <- methylEst(bm, controlCI = list(compute = TRUE, method = "Wald", 
            level = 0.95, nmarg = 512, ncpu = NULL))
    }
}
\keyword{programming}
