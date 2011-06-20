\name{enrichmentCalc}
\alias{enrichmentCalc}
\alias{enrichmentCalc,GRanges-method}
\alias{enrichmentCalc,GRangesList-method}
\title{Calculate sequencing enrichment}
\description{Function to calculate enrichment over the whole genome of sequencing reads.}

\usage{
  \S4method{enrichmentCalc}{GRanges}(x, seq.len = NULL, verbose = TRUE)
  \S4method{enrichmentCalc}{GRangesList}(x, verbose = TRUE, ...)
}
\arguments{
  \item{x}{A \code{GRangesList} or \code{GRanges} object. All chromosome lengths
           must be stored in the \code{Seqinfo} of this object.}
  \item{seq.len}{If sequencing reads need to be extended, the fragment size to be used.}
  \item{verbose}{Whether to print the progress of processing.}
  \item{...}{Argument \code{seq.len} above, not directly used in the \code{GRangesList} method.}
}

\details{
  If \code{seq.len} is supplied, \code{x} is firstly extended, and then turned into
  a coverage object. The number of extended reads covering each base pair of the
  genome is then tabulated, and returned as a \code{data.frame}.
}
\value{
  For the \code{GRanges} method, \code{data.frame} containing columns
  \code{coverage} and \code{bases}. For the \code{GRangesList} method,
  a list of such \code{data.frame}s.
}
\author{Aaron Statham}
\examples{
  require(GenomicRanges)
  readsFile <- system.file("data", "samplesList.RData", package = "Repitools")
  load(readsFile)  # GRangesList of reads 'samplesList'
  seqlengths(samplesList)

  tc <- enrichmentCalc(samplesList, seq.len = 300)
}
