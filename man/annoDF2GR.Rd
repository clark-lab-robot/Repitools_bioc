\name{annoDF2GR}
\alias{annoDF2GR}
\alias{annoDF2GR,data.frame-method}
\title{Convert a \code{data.frame} to a \code{GRanges}.}
\description{
  Checks that the \code{data.frame} has the required columns, \code{chr},
  \code{start}, \code{end}, then creates a \code{GRanges}, keeping all of the
  additional columns.
}
\usage{
  \S4method{annoDF2GR}{data.frame}(anno)
}
\arguments{
  \item{anno}{An \code{data.frame}, describing some genomic features.}
}
\details{
  Extra columns are added to the \code{elementMetadata} of the \code{GRanges} object.
}
\value{
  A \code{\linkS4class{GRanges}} of the annotation.
}
\author{Dario Strbenac}
\examples{
  df <- data.frame(chr = c("chr1", "chr3", "chr7", "chr22"),
                   start = seq(1000, 4000, 1000),
                   end = seq(1500, 4500, 1000),
                   t = c(3.11, 0.93, 2.28, -0.18),
                   gc = c("High", "High", "Low", "High"))

  annoDF2GR(df)
}
