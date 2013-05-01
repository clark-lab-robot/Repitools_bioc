\name{empBayes}
\alias{empBayes}
\title{
    Function to calculate prior parameters using empirical Bayes.
}
\description{
    Under the empirical Bayes approach (and assuming a uniform prior for 
    the methylation level) the shape and scale parameters for the gamma prior
    of the region-specific read density are derived.  The parameters are 
    thereby determined in a CpG-dependent manner.
}
\usage{
empBayes(x, ngroups = 100, ncomp = 1, maxBins=50000, method="beta", controlMethod=list(mode="full", weights=c(0.1, 0.8, 0.1), param=c(1,1)), ncpu = NULL, verbose = FALSE)
}
\arguments{
    \item{x}{
    Object of class \code{BayMethList}.
}
    \item{ngroups}{
    Number of CpG density groups you would like to consider. The bins are
    classified based on its CpG density into one of \code{ngroups} classes and
    for each class separately the set of prior parameters will be determined.
}
    \item{ncomp}{
    Number of components of beta distributions in the prior distribution for
    the methylation level when method is equal to \code{beta}. 
}
    \item{maxBins}{
    Maximum number of bins in one CpG density group used to derive the 
    parameter estimates. If maxBins is smaller than the number of bins
    that are in one groups than \code{maxBins} bins are sampled with 
    replacement.
}
    \item{method}{
    Either \code{DBD} for a Dirac-Beta-Dirac mixture, representing a mixture
    a mixture of a point mass at zero, a beta distribution and a point mass at one,
    or \code{beta} for as Beta mixture with \code{ncomp} components.
}
  \item{controlMethod}{
    list defining settings if the Dirac-Beta-Dirac mixture is chosen. 
    \describe{
    \item{-}{
        \code{mode}
        Either \code{full}, \code{fixedWeights} or \code{fixedBeta}. 
        Using the \code{full} both the mixture weights and 
        beta parameters are estimated. In mode \code{fixedWeights} the weights
        are fixed given to the values in \code{weights} and only the parameters
        of the beta component are estimated. 
        In mode \code{fixedBeta} the parameters of the beta component are 
        fixed to the values specified in \code{param}. 
        The default mode is \code{full}.
    }
    \item{-}{
        \code{weights}
        Numeric vector of length three specifying the weights for the Dirac-Beta-Dirac
        mixture when mode is equal to \code{fixedWeights}. 
        The first element specifies the weight for the zero point mass, the
        second for the beta component and the third for the point mass at one.
        The three values must sum up to one. The default is c(0.1, 0.8, 0.1).
    }
    \item{-}{
        \code{param}
        Numeric vector of length two specifying (positive) parameters of the beta distribution 
        component when mode is equal to \code{fixedBeta}. 
        The default is c(1,1).
    }
  }
}
    \item{ncpu}{
    Number of CPUs on your machine you would like to use in parallel. 
    If \code{ncpu} is set to NULL, half of the CPUs will be used on machines
    with a maximum of four CPUs, and 2/3 will be used if more CPUs are available.
}
    \item{verbose}{
    Boolean indicating whether the empirical Bayes function should run in a verbose mode (default 'FALSE').
} 
}
\details{
BayMeth takes advantage of the relationship between CpG-density and read 
depth to formulate a CpG-density-dependent gamma prior distribution for the
region-specific read density. Taking CpG-density into account the prior should
stabilise the methylation estimation procedure for low counts and in the
presence of sampling variability. The shape and scale parameter of the gamma
prior distribution are determined in a CpG-density-dependent manner using
empirical Bayes. For each genomic bin the CpG density is provided in the
\code{BayMethList}-object. Each bin is classified based on its CpG-density into
one of \code{ngroups} non-overlapping CpG-density intervals. For each class
separately, we derive the values for the shape and scale parameter under an
empirical Bayes framework using maximum likelihood. For CpG classes which
contain more than \code{maxBins} bins, a random sample drawn with replacement of size
\code{maxBins} is used to derive these prior parameters. Note that both read depths,
from the SssI control and the sample of interest, are thereby taken into
account. We end up with \code{ngroups} parameter sets for shape and rate.
}
\value{
A \code{BayMethList} object where the slot \code{priorTab} is filled. \code{priorTab} represent a list. The first list entry contains the CpG group a bin is assigned to. The second entry contains the number of components that have been used for the prior (at the moment 1).  The following list entries correspond to one sample of interest, respectively, and contain a matrix with the optimal shape and scale parameters for all CpG classes. The first row contains the optimal shape parameter and the second row the optimal scale parameter. The number of columns corresponds to the number of CpG classes specified in \code{ngroups}.
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
 
        # mask out unannotated high copy number regions
        # see Pickrell et al. (2011), Bioinformatics 27: 2144-2146.

        # should take about 3 minutes for both sample of interests with 2 CPUs.
        bm <- empBayes(bm, ngroups=20) 
   }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{programming}
