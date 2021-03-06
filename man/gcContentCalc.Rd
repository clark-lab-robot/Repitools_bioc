\name{gcContentCalc}
\alias{gcContentCalc}
\alias{gcContentCalc,GRanges,BSgenome-method}
\alias{gcContentCalc,data.frame,BSgenome-method}
\title{Calculate The gcContent of a Region}
\description{Function to calculate the GC content of windows}
\usage{
  \S4method{gcContentCalc}{GRanges,BSgenome}(x, organism, verbose = TRUE)
  \S4method{gcContentCalc}{data.frame,BSgenome}(x, organism, window = NULL, ...)
}
\arguments{
  \item{x}{A \code{GRanges} object or a \code{data.frame}, with columns \code{chr}
           and either \code{position} or \code{start}, \code{end} and \code{strand}.}
  \item{window}{Bases around the locations that are in the window. Calculation will
                consider \code{windowSize/2} bases upstream, and \code{windowSize / 2 - 1 }
                bases downstream.}
  \item{organism}{The \code{BSgenome} object to calculate gcContent upon.}
  \item{verbose}{Whether to print the progess of processing.}
  \item{...}{The \code{verbose} variable for the \code{data.frame} method,
             passed onto the \code{GRanges} method.}
}
\details{
  The windows considered will be \code{windowSize/2} bases upstream and \code{windowSize/2-1}
  bases downstream of the given position, for each position. The value returned for
  each region is a percentage of bases in that region that are a G or C.
}
\value{
  A vector of GC content percentages, one for each region.
}
\author{Aaron Statham}
\examples{
require(BSgenome.Hsapiens.UCSC.hg18)
TSSTable <- data.frame(chr = paste("chr", c(1,2), sep = ""), position = c(100000, 200000))
gcContentCalc(TSSTable, 200, organism=Hsapiens)
}

