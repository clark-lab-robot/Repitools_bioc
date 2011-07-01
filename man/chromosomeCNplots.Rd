\name{chromosomeCNplots}
\alias{chromosomeCNplots}
\alias{chromosomeCNplots,CopyEstimate-method}
\title{Plot GC content vs. Read Counts Before Normalising, and GC content vs. Copy Estimates
       After Normalising.}
\description{
  Two plots on the same plotting page are made for each sample. The top plot has estimates
  of copy number for a chromosome before any GC correction was applied. The bottom
  plot shows the copy number estimates for a chromosome after GC correction was applied.
}
\usage{
  \S4method{chromosomeCNplots}{CopyEstimate}(copy, y.max = NULL, pch = 19, cex = 0.2,
             pch.col = "black", seg.col = "red", lty = 1, lwd = 2, verbose = TRUE)
}
\arguments{
  \item{copy}{A \code{\linkS4class{CopyEstimate}} object.}
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
  copy number estimates. The segmentation line plotted in the bottom figure on each page
  is the segmentation regions found.
}
\value{
  A number of pages of scatterplots. Each page has two plots for a chromosome,
  one for before GC correction and one for after. The output should, therefore, be sent
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


