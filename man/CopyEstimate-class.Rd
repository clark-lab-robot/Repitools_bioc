\name{CopyEstimate}
\docType{class}
\alias{CopyEstimate}
\alias{CopyEstimate-class}
\alias{CopyEstimate,GRanges,matrix,GRangesList,character-method}
\alias{show,CopyEstimate-method}

\title{Container for results of fold change copy number estimation.}
\description{
  Contains the genomic coordinates of regions, and fold change estimates.
}
\section{Constructor}{
  \describe{
    \item{}{
    	\code{CopyEstimate(windows, unadj.CN, unadj.CN.seg)}
	Creates a CopyEstimate object.
	\describe{
	    \item{\code{windows}}{A \code{\linkS4class{GRanges}} object.}
	    \item{\code{unadj.CN}}{A matrix of fold changes.}
	    \item{\code{unadj.CN.seg}}{A \code{\linkS4class{GRangesList}} object holding
                                       the segmentation results.}
	}
    }
  }
}
\section{Additional Slots}{
  \describe{
      These are added to by \code{\link{absoluteCN}} or \code{\link{relativeCN}}
    \item{type}{A flag that contains if the copy number data is absolute or relative.}
  }
}

