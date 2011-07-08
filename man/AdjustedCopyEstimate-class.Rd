\name{AdjustedCopyEstimate}
\docType{class}
\alias{AdjustedCopyEstimate}
\alias{AdjustedCopyEstimate-class}
\alias{AdjustedCopyEstimate,numeric,GRanges,numeric,numeric,matrix,list,matrix,character-method}
\alias{show,AdjustedCopyEstimate-method}

\title{Container for results of GC adjusted copy number estimation.}
\description{
  Contains the genomic coordinates of regions, the raw counts before GC adjustment,
  the GC content and mappability of each region, and the polynomial model fit, and
  the GC-adjusted copy number estimates.
}
\section{Constructor}{
  \describe{
    \item{}{
    	\code{AdjustedCopyEstimate(ploidy, windows, mappability, gc, unadj.CN, models, adj.CN)}
	Creates a AdjustedCopyEstimate object.
	\describe{
            \item{\code{ploidy}}{Sets of chromosomes in each sample.}
	    \item{\code{windows}}{A \code{\linkS4class{GRanges}} object.}
	    \item{\code{mappability}}{A numeric vector of mappability. Elements
                                      between 0 and 1.}
	    \item{\code{gc}}{A numeric vector of GC content Elements between 0 and 1.}
	    \item{\code{unadj.CN}}{A matrix of estimated copy numbers after mappability
                                   adjustment, but before GC content adjustment, if slot \code{type}
                                   is \code{"absolute"}. Otherwise, fold changes.}
	    \item{\code{models}}{The polynomial models that were fit to the counts.}
	    \item{\code{adj.CN}}{A matrix of estimated copy numbers after mappability
                                   adjustment and GC content adjustment, if slot \code{type}
                                   is \code{"absolute"}. Otherwise, a matrix of fold changes,
                                   based on GC adjusted absolute copy estimates.}
	}
    }

    Note that \code{mappability} and \code{gc} become metadata columns of \code{windows}
    when the object is created.
  }
}

\section{Superclass}{
  \describe{
        This class inherits from \code{\linkS4class{CopyEstimate}}.
    }
}

\section{Additional Slots}{
  \describe{
      These are added to by \code{\link{absoluteCN}} or \code{\link{relativeCN}}
    \item{unadj.CN.seg}{A \code{\linkS4class{GRangesList}} of copy number segmentations
                  for each sample.}
    \item{adj.CN.seg}{A \code{\linkS4class{GRangesList}} of copy number segmentations
                  for each sample, using GC adjusted data.}
    \item{type}{A flag that contains if the copy number data is absolute or relative.}
  }
}
