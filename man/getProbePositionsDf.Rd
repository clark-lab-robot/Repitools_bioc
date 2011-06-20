\name{getProbePositionsDf}
\alias{getProbePositionsDf}
\alias{getProbePositionsDf,AffymetrixCdfFile-method}
\title{Translate Affymetrix probe information in a table.}
\description{
  Translates the probe information in the AromaCellPositionFile
  to a data.frame object.
}
\usage{
  \S4method{getProbePositionsDf}{AffymetrixCdfFile}(cdf, chrs, ..., verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cdf}{An AffymetrixCdfFile object.}
  \item{chrs}{A vector of chromosome names. Optional.}
  \item{\dots}{Further arguments to send to \code{getCellIndices}.}
  \item{verbose}{Logical; whether or not to print out progress statements to the screen.}
}
\details{
This assumes that the AromaCellPositionFile exist.
}
\value{
  A data.frame with 3 columns: chr, position, index
}
\author{Mark Robinson}
\examples{
## not run
# probePositions <- getProbePositionsDf(cdfU)
}
