\name{absoluteCN}
\alias{absoluteCN}
\alias{absoluteCN,data.frame,matrix,GCAdjustParams-method}
\alias{absoluteCN,GRanges,matrix,GCAdjustParams-method}
\title{Calculate and Segment Absolute Copy Number from Sequencing Counts}
\description{This function uses the \code{\link{GCadjustCopy}} function to convert
  a matrix of count data into absolute copy number estimates, then it segments them,
  and reports the copy number of either the input regions or user-defined regions of
  interest.
}
\usage{
  \S4method{absoluteCN}{data.frame,matrix,GCAdjustParams}(input.windows, input.counts, gc.params, ...)
  \S4method{absoluteCN}{GRanges,matrix,GCAdjustParams}(input.windows, input.counts, gc.params,
                                                       segment.sqrt = TRUE, ..., verbose = TRUE)
}
\arguments{
  \item{input.windows}{A \code{data.frame} with (at least) columns \code{chr},
                       \code{start}, and \code{end}, or a GRanges object.}
  \item{input.counts}{A matrix of counts. Rows are genomic windows and columns are
                      samples.}
  \item{gc.params}{A \code{\linkS4class{GCAdjustParams}} object, holding parameters
                   related to mappability and GC content correction of read counts.}
  \item{segment.sqrt}{Whether to square root the absolute copy number estimates before
                      running the segmentation.}
  \item{...}{For the \code{data.frame} method; the \code{verbose} variable and any
             additional parameters to pass to the \code{segment} function. For the
             \code{GRanges} method; additional parameters for the segmentation.}
  \item{verbose}{Whether to print the progess of processing.}
}
\details{
  For details of the absolute copy number estimation step, see the documentation for
  \code{\link{GCadjustCopy}}.

  For details of the segmentation, see \code{\link[DNAcopy]{segment}} documentation.
  By default, no weights are used.
}
\value{
  A \code{\linkS4class{CopyEstimate}} object. If \code{regions} was not provided,
  it describes the input windows, otherwise it describes the windows specified by
  \code{regions}.
}
\author{Dario Strbenac}
\examples{
  \dontrun{
    library(BSgenome.Hsapiens.UCSC.hg18)
    library(BSgenome.Hsapiens36bp.UCSC.hg18mappability)
    load("inputsReads.RData")
    windows <- genomeBlocks(Hsapiens, chrs = paste("chr", c(1:22, 'X', 'Y'), sep = ''),
                            width = 20000)
    counts <- annotationBlocksCounts(inputsReads, anno = windows, seq.len = 300)

    gc.par <- GCAdjustParams(genome = Hsapiens, mappability = Hsapiens36bp,
                             min.mappability = 50, n.bins = 10, min.bin.size = 10,
                             poly.degree = 4, ploidy = c(2, 4))
    abs.cn <- absoluteCN(input.windows = windows, input.counts = counts, gc.params = gc.par)
  }
}

