\name{annoGR2DF}
\alias{annoGR2DF}
\alias{annoGR2DF,GRanges-method}
\title{Convert an annotated \code{GRanges} to a \code{data.frame}.}
\description{
  Converting a \code{GRanges} that might be annotated with some kind of results
  to a \code{data.frame} is useful, because it allows easier writing to file
  and viewing in other programs, like a spreadsheet program.
}
\usage{
  \S4method{annoGR2DF}{GRanges}(anno)
}
\arguments{
  \item{anno}{A \code{GRanges}, describing some genomic features.}
}
\details{
  The column name \code{seqnames} is changed to \code{chr}, and if all the
  strands are \code{*}, then the \code{strand} column is dropped.
}
\value{
  A \code{data.frame} of the annotation.
}
\author{Dario Strbenac}
\examples{
  require(GenomicRanges)
  chrs <- c("chr1", "chr3", "chr7", "chr22")
  starts <- seq(1000, 4000, 1000)
  ends <- seq(1500, 4500, 1000)
  t <- c(3.11, 0.93, 2.28, -0.18)
  gc <- c("High", "High", "Low", "High")
  gr <- GRanges(chrs, IRanges(starts, ends), strand = '*', t, gc)

  annoGR2DF(gr)
}
