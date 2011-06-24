\name{absoluteCN}
\alias{absoluteCN}
\alias{absoluteCN,data.frame,data.frame-method}
\alias{absoluteCN,GRanges,GRanges-method}
\title{Calculate Absolute Copy Number from Sequencing Counts}
\description{Taking into account mappability and GC content biases, the absolute
  copy number is calculated, by assuming that the median read depth is a copy number
  of 1.
}
\usage{
  \S4method{absoluteCN}{data.frame,data.frame}(input.windows, ip.windows, ...)
  \S4method{absoluteCN}{GRanges,GRanges}(input.windows, ip.windows, input.counts = NULL,
                                         gc.params = NULL, verbose = TRUE)
}
\arguments{
  \item{input.windows}{A \code{data.frame} with (at least) columns \code{chr},
                       \code{start}, and \code{end}, or a GRanges object.}
  \item{ip.windows}{A \code{data.frame} with (at least) columns \code{chr},
                       \code{start}, and \code{end}, or a GRanges object.}
  \item{input.counts}{A matrix of counts. Rows are genomic windows and columns are
                      samples.}
  \item{gc.params}{A \code{\linkS4class{GCAdjustParams}} object, holding parameters
                   related to mappability and GC content correction of read counts.}
  \item{...}{Parameters not directly used by the conversion of the \code{data.frame}
             to a \code{GRanges} object.}
  \item{verbose}{Whether to print the progess of processing.}
}
\details{
  First, the mappability of all counting windows is calculated, and windows that
  have mappability less than the cutoff specified by in the parameters object are
  ignored in further steps. The remaining windows have their counts scaled by
  multiplying their counts by 100 / percentage mappability.

  The range of GC content of the counting windows is broken into a number of bins,
  as specified by the user in the parameters object. A probability density function
  is fitted to the counts in each bin, so the mode can be found. The mode is taken
  to be the counts of the copy neutral windows, for that GC content bins.

  A polynomial function is fitted to the modes of GC content bins. Each count is
  divided by its expected counts from the polynomial function to give an absolute
  copy number estimate. If the ploidy has been provided in the parameters object,
  then all counts within a sample are multiplied by the ploidy for that sample.
  If the sample ploidys were omitted, then no scaling for ploidy is done.
}
\value{
  A matrix of copy numbers of the regions specified by \code{ip.windows}. Each column
  of the matrix is one sample. If region of interest did not have an estimated copy
  number, it will contain NA.
}
\author{Dario Strbenac}
\examples{
  \dontrun{
    library(BSgenome.Hsapiens.UCSC.hg18)
    library(BSgenome.Hsapiens36bp.UCSC.hg18mappability)
    load("inputsReads.RData")
    windows <- genomeBlocks(Hsapiens, chrs = paste("chr", c(1:22, 'X', 'Y'), sep = ''),
                            width = 20000)
    counts <- annotationBlocksCounts(inputsReads, windows, seq.len = 300)

    gc.par <- GCAdjustParams(Hsapiens, Hsapiens36bp, 50, 10, 10, 4, c(2, 4))
    abs.cn <- absoluteCN(windows, windows, counts, gc.par)
  }
}

