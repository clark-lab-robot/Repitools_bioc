\name{CopyEstimate}
\docType{class}
\alias{CopyEstimate}
\alias{CopyEstimate-class}
\alias{CopyEstimate,GRanges,matrix,numeric,numeric,list,matrix-method}
\alias{show,CopyEstimate-method}

\title{Container for results of GC adjusted copy number estimation.}
\description{
  Contains the genomic coordinates of regions, the raw counts before GC adjustment,
  the GC content and mappability of each region, and the polynomial model fit, and
  the GC-adjusted copy number estimates.
}
\section{Constructor}{
  \describe{
    \item{}{
    	\code{CopyEstimate(windows, old.counts, mappability, gc, models, cn)}
	Creates a CopyEstimate object.
	\describe{
	    \item{\code{windows}}{A \code{\linkS4class{GRanges}} object.}
	    \item{\code{old.counts}}{A matrix of counts before any adjustment.}
	    \item{\code{mappability}}{A numeric vector of mappability. Elements
                                      between 0 and 1.}
	    \item{\code{gc}}{A numeric vector of GC content Elements between 0 and 1.}
	    \item{\code{models}}{The polynomial model that was fit to the counts.}
	    \item{\code{cn}}{A matrix of GC-adjusted absolute copy number estimates.}
	}
    }
  }
}
