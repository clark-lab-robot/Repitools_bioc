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
methylEst(x, verbose=FALSE, controlCI = list(compute = FALSE, method = "Wald", 
    level = 0.95, nmarg = 512, ncpu = NULL))
}
\arguments{
  \item{x}{
    Object of class \code{BayMethList}.
    }
    \item{verbose}{
    Boolean indicating whether the methylEst function should run in a verbose mode (default 'FALSE').
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
        There are three possible types of credible intervals that can be chosen if a uniform prior, e.g. Beta(1,1), is chosen:
        'Wald' (default), 'HPD', 'quantile'. The Wald-type intervals are the fastest to compute. The are calculated on logit scale 
        and then transferred back. Due to numerical integration of the posterior marginal posterior distributions, the computation of highest posterior density (HPD) interval and quantile-based interval is computationally more expensive. However, in our applications HPD intervals provided best coverage. 

        Note, using a beta mixture or a Dirac-beta-Dirac (DBD) mixture as prior distribution for the methylation level only \code{method="quantile"} is available.
    }
    \item{-}{
        \code{level}
        numerical value defining the credible level. Default: 0.95.
        }
    \item{-}{
        \code{nmarg}
        Number of points at which the posterior marginal is evaluated
        (only relevant for method="quantile" or method="HPD").
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
straightforward to efficiently compute; Wald-based credible intervals are obtained
on logit scale and then back-transferred to ensure values withing 0 and 1. HPD and
quantile-based credible intervals are computed by numerical integration of the posterior marginal distribution.
}
\value{
    A \code{BayMethList} object where the slot \code{methEst} is filled with a
    list containing the following elements:
    \item{mean }{Matrix where the number of columns equals the number of
    samples of interest. Each column contains the posterior mean methylation
    level for each bin.}
    \item{var }{Matrix where the number of columns equals the number of samples 
        of interest. Each column contains posterior variance for each bin.}
    \item{ci }{List with length equal to the number of samples of interest. Each list element contains a matrix where the first column contains the lower CI bound, the second column the upper CI bound, and the third column the CI width.}
    \item{W }{Matrix where the number of columns equals the number of
    samples of interest. Each column contains the normalisation factor of the posterior marginal distribution 
        for each bin.}
    \item{al }{Matrix where the number of columns equals the number of
    samples of interest. Each column contains the prior shape parameter for each bin}
    \item{bl }{Matrix where the number of columns equals the number of
    samples of interest. Each column contains the prior scale parameter for each bin}
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
        bm <- methylEst(bm, controlCI = list(compute = FALSE, method = "Wald", 
            level = 0.95, nmarg = 512, ncpu = NULL))
    }
}
\keyword{programming}
