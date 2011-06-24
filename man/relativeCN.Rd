\name{relativeCN}
\alias{relativeCN}
\alias{relativeCN,data.frame,data.frame-method}
\alias{relativeCN,GRanges,GRanges-method}
\title{Calculate Relative Copy Number and Map To Enriched Regions}
\description{Based on the regions defined in \code{input.windows}, the function
  calculates the log2 fold change between the treatment and control columns
  of counts. Then, based on the regions defined in the \code{ip.windows} element,
  it finds which copy number region each enrichment region belongs to, and returns
  a vector that holds the copy number state of every enrichment region.
}
\usage{
  \S4method{relativeCN}{data.frame,data.frame}(input.windows, ip.windows, input.counts = NULL,
                                               gc.params = NULL, ..., verbose = TRUE)
  \S4method{relativeCN}{GRanges,GRanges}(input.windows, ip.windows, input.counts = NULL,
                                         gc.params = NULL, ..., verbose = TRUE)
}
\arguments{
  \item{input.windows}{A \code{data.frame} with (at least) columns \code{chr},
                       \code{start}, and \code{end}, or a GRanges object.}
  \item{ip.windows}{A \code{data.frame} with (at least) columns \code{chr},
                    \code{start}, and \code{end}, or a GRanges object.}
  \item{input.counts}{A matrix of counts. The first column must be for the control
                      state, and the second column must be for the treatment state.}
  \item{gc.params}{A \code{\linkS4class{GCAdjustParams}} object, holding parameters
                   related to mappability and GC content correction of read counts.}
  \item{...}{Further parameters passed to \code{\link{segment}} function in \code{DNAcopy}
             package.}
  \item{verbose}{Whether to print the progess of processing.}
}
\details{
  The algorithm used to decide on the copy number regions is Circular Binary
  Segmentation (Olshen et al. 2004). If an enrichment region maps to a region
  that was not segmentable, it is assigned NA.

  If \code{gc.params} is \code{NULL}, then no correction for mappability
  or GC content is done. If it is specified, then absolute copy numbers are
  estimated with \code{\link{absoluteCN}} for each condition, which corrects
  for mappability and then GC content, before estimating absolute copy numbers.
  The ratio of estimated absolute copy numbers is used to calculate relative
  copy numbers.
}
\value{
  A vector with the same number of elements as \code{input.windows}, that has
  the copy number of each region relative to the control. The copy number ratios
  returned are on the linear scale, not log2.
}
\author{Dario Strbenac}
\references{Olshen, A. B., Venkatraman, E. S., Lucito, R., and Wigler, M. (2004).
  Circular binary segmentation for the analysis of array-based DNA copy number data.
  \emph{Biostatistics} 5: 557-572}
\examples{
  inputs <- data.frame(chr = c("chr1", "chr1", "chr1", "chr2", "chr2"),
                     start = c(1, 50001, 100001, 1, 10001),
                       end = c(50000, 100000, 150000, 10000, 20000))
  enriched <- data.frame(chr = c("chr1", "chr1", "chr1", "chr1", "chr1", "chr2", "chr2"),
                       start = c(1, 1001, 2001, 10000001, 10005001, 9001, 9501),
                         end = c(1000, 2000, 3000, 10001000, 10006000, 9500, 10000))
  counts <- matrix(c(25, 39, 3, 10, 22, 29, 38, 5, 19, 31), nrow = 5)
  colnames(counts) <- c("Control", "Treatment")
  relativeCN(inputs, enriched, counts, p.method = "perm")
}

