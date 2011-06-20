\name{enrichmentPlot}
\alias{enrichmentPlot}
\alias{enrichmentPlot,GRangesList-method}
\title{Plot the distribution of sequencing enrichment.}
\description{Function to generate a plot of the distribution of sequencing reads enrichments.}

\usage{
  \S4method{enrichmentPlot}{GRangesList}(x, seq.len, cols = rainbow(length(x)),
      xlim = c(0, 20), main = "Enrichment Plot", total.lib.size = TRUE, verbose = TRUE, ...)
}

\arguments{
  \item{x}{A \code{GRangesList} object of reads to plot enrichment of. The chromosome
           lengths must be stored in the \code{Seqinfo} of this object.}
  \item{seq.len}{The fragment size to be used for extending the sequencing reads.}
  \item{cols}{The line colour for each element of \code{x}}
  \item{xlim}{\code{xlim} parameter passed to \code{plot}, the default is appropriate for
              \code{"linear"} \code{cpgDensityCalc} weighting.}
  \item{main}{\code{main} parameter passed to \code{plot}}
  \item{total.lib.size}{Whether to normalise enrichment values to the total number
                        of reads per lane.}
  \item{verbose}{Print details of processing.}
  \item{...}{Additional graphical parameters to pass to \code{plot}.}
}

\details{
    See \code{enrichmentCalc} for details of how the results are determined.
}

\value{
    A plot is created. The data processed by \code{enrichmentCalc} is invisibly returned.
}

\author{Aaron Statham}
\examples{
  readsFile <- system.file("data", "samplesList.RData", package = "Repitools")
  load(readsFile)  # GRangesList of reads 'samplesList'
  enrichmentPlot(samplesList, seq.len = 300, total.lib.size = FALSE)
}
