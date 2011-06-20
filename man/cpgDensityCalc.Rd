\name{cpgDensityCalc}
\alias{cpgDensityCalc}
\alias{cpgDensityCalc,data.frame,BSgenome-method}
\alias{cpgDensityCalc,GRangesList,BSgenome-method}
\alias{cpgDensityCalc,GRanges,BSgenome-method}
\title{Calculate CpG Density in a Window}
\description{Function to calculate CpG density around a position.
}
\usage{
  \S4method{cpgDensityCalc}{data.frame,BSgenome}(x, organism, ...)
  \S4method{cpgDensityCalc}{GRangesList,BSgenome}(x, organism, verbose = TRUE, ...)
  \S4method{cpgDensityCalc}{GRanges,BSgenome}(x, organism, seq.len = NULL, window = NULL,
                                          w.function = c("none", "linear", "exp", "log"),
                                          verbose = TRUE)
}
\arguments{
  \item{x}{A \code{data.frame}, with columns \code{chr} and \code{position},
           or columns \code{chr}, \code{start}, \code{end}, and \code{strand}.
           Also may be a \code{GRangesList} object, or \code{GRanges}.}
  \item{window}{Bases around the locations that are in the window. Calculation
               will consider \code{window/2 - 1} bases upstream, and
               \code{window/2} bases downstream.}
  \item{w.function}{Weighting function to use. Can be \code{"none"}, \code{"linear"},
                    \code{"log"}, or \code{"exp"}}
  \item{organism}{The \code{BSgenome} object to calculate CpG density upon.}
  \item{seq.len}{The fragment size of the sequence reads in \code{x}.
                 Default: No extension.}
  \item{verbose}{Print details of processing.}
  \item{...}{Arguments passed into the \code{data.frame} or \code{GRangesList} method,
             but not used until the \code{GRanges} method.}
}
\details{
  If the version of the data frame with the start, end, and strand columns is given,
  the window will be created around the TSS.
  
  For weighting scheme \code{"none"}, this is equivalent to the number of CG matches
  in the region. For \code{"linear"} weighting, each match is given a score
  \code{1/x} where \code{x} is the number of bases from the postition that the match occurred,
  and the scores are summed. For exponential weighting and logarithmic weighting,
  the idea is similar, but the scores decay exponentially (\code{exp^-5x/window})
  and logarithmically (\code{log2(2 - (distancesForRegion / window))}).
}
\value{
  A \code{numeric} vector of CpG densities for each region.
}
\author{Dario Strbenac}
\examples{
  if(require(BSgenome.Hsapiens.UCSC.hg18))
  {
    TSSTable <- data.frame(chr = c("chr1", "chr2"), position = c(100000, 200000))
    cpgDensityCalc(TSSTable, organism = Hsapiens, window = 600)
  }
}
