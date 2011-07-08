\name{chromosomeCNplots}
\alias{chromosomeCNplots}
\alias{chromosomeCNplots,CopyEstimate-method}
\alias{chromosomeCNplots,AdjustedCopyEstimate-method}
\title{Plot copy number by chromosome}
\description{
  Generates plots of position along chromosomes vs. estimated copy number. If GC adjustment was
  performed, then there are two plots per page; one before adjustment and one after adjustment.
}
\usage{
  \S4method{chromosomeCNplots}{CopyEstimate}(copy, y.max = NULL, pch = 19, cex = 0.2,
             pch.col = "black", seg.col = "red", lty = 1, lwd = 2, verbose = TRUE)
  \S4method{chromosomeCNplots}{AdjustedCopyEstimate}(copy, y.max = NULL, pch = 19, cex = 0.2,
             pch.col = "black", seg.col = "red", lty = 1, lwd = 2, verbose = TRUE)
}
\arguments{
  \item{copy}{A \code{\linkS4class{CopyEstimate}} or \code{\linkS4class{AdjustedCopyEstimate}} object.}
  \item{y.max}{The maximum value of the y-axis of the scatter plots.}
  \item{pch}{Style of points in the scatter plots.}
  \item{cex}{Whether to square root the absolute copy number estimates before
                      running the segmentation.}
  \item{pch.col}{Colour of points in the scatter plots.}
  \item{seg.col}{Colour of copy number segmentation line.}
  \item{lty}{Line type of plotted regression line.}
  \item{lwd}{Line width of plotted regression line.}
  \item{verbose}{Whether to print the progess of processing.}
}
\details{
  See \code{\link{absoluteCN}} or \code{\link{relativeCN}} for how to do the GC adjusted
  copy number estimates, if this is required. The segmentation line plotted
  is of the segmentation regions found by circular binary segmentation.
}
\value{
  A number of pages of scatterplots. The output should, therefore, be sent
  to a PDF device.
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

    pdf("chrProfiles.pdf")
    chromosomeCNplots(abs.cn, y.max = 8)
    dev.off()
  }
}


