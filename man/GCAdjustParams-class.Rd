\name{GCAdjustParams}
\docType{class}
\alias{GCAdjustParams}
\alias{GCAdjustParams-class}
\alias{GCAdjustParams,BSgenome,BSgenome-method}

\title{Container for parameters for mappability and GC content adjusted absolute
       copy number estimation.}
\description{
  The parameters are used by the \code{\link{absoluteCN}} function.
}
\section{Constructor}{
  \describe{
    \item{}{
    	\code{GCAdjustParams(genome, mappability, min.mappability, n.bins = NULL,
                             min.bin.size = 1, poly.degree = NULL, ploidy = 1)}
	Creates a GCAdjustParams object.
	\describe{
	    \item{\code{genome}}{A \code{\linkS4class{BSgenome}} object of the
                                 species that the experiment was done for.}
	    \item{\code{mappability}}{A \code{\linkS4class{BSgenome}} object, containing
                                      the mappability of each base in the genome.}
	    \item{\code{min.mappability}}{A number between 0 and 100 that is a cutoff
                                          on window mappability.}
	    \item{\code{n.bins}}{The number of GC content bins to divide the windows
                                 into, before finding the mode of counts in each window.}
	    \item{\code{min.bin.size}}{GC bins with less than this many count windows inside
                                       them will be ignored.}
	    \item{\code{poly.degree}}{The degree of the polynomial to fit to the GC bins'
                                      count modes.}
	    \item{\code{ploidy}}{A vector of multipliers to use on the estimated
                                absolute copy number of each sample,
                                if the number of sets of chromosomes is known.}
	}
    }
  }
}
\author{Dario Strbenac}

