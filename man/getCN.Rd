\name{getCN}
\alias{getCN}
\alias{getCN,data.frame,data.frame-method}
\alias{getCN,GRanges,GRanges-method}
\title{Calculate Copy Number and Map To Enriched Regions}
\description{Based on the regions defined in \code{input.windows}, the function
  first calculates the log2 fold change between the treatment and control columns
  of counts. Then, based on the regions defined in the \code{ip.windows} element,
  it finds which copy number region each enrichment region belongs to, and returns
  a vector that holds the copy number state of every enrichment region.
}
\usage{
  \S4method{getCN}{data.frame,data.frame}(input.windows, ip.windows, ...)
  \S4method{getCN}{GRanges,GRanges}(input.windows, ip.windows, input.counts = NULL,
                                    ..., verbose = TRUE)
}
\arguments{
  \item{input.windows}{A \code{data.frame} with (at least) columns \code{chr},
                       \code{start}, and \code{end}, or a GRanges object.}
  \item{ip.windows}{A \code{data.frame} with (at least) columns \code{chr},
                    \code{start}, and \code{end}, or a GRanges object.}
  \item{input.counts}{A matrix of counts. The first column must be for the control
                      state, and the second column must be for the treatment state.}
  \item{...}{Further parameters passed to getCN. If the \code{data.frame,data.frame}
             method is called, \code{input.counts} and \code{verbose} will also be
             passed in this way.}
  \item{verbose}{Whether to print the progess of processing. Default: TRUE.}
}
\details{
  The algorithm used to decide on the copy number regions is Circular Binary
  Segmentation (Olshen et al. 2004). If an enrichment region maps to a region
  that was not found to be in a copy number region from the input samples,
  it is assigned the copy number of 1.
}
\value{
  A vector with the same number of elements as \code{input.windows}, that has
  the copy number of each region relative to the control.
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
  getCN(inputs, enriched, counts, p.method = "perm")
}

