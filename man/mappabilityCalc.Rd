\name{mappabilityCalc}
\alias{mappabilityCalc}
\alias{mappabilityCalc,GRanges-method}
\alias{mappabilityCalc,data.frame-method}
\title{Calculate The Mappability of a Region}
\description{Function to calculate mappability of windows}
\usage{
  \S4method{mappabilityCalc}{GRanges}(x, organism, window = NULL,
              type = c("TSS", "center", "block"), verbose = TRUE)
  \S4method{mappabilityCalc}{data.frame}(x, organism, window = NULL,
                                         type = c("TSS", "center", "block"), ...)
}
\arguments{
  \item{x}{A \code{GRanges} object or a \code{data.frame}, with columns \code{chr} and
           either \code{position} or \code{start}, \code{end} and \code{strand}.}
  \item{window}{Bases around the locations that are in the window. Calculation will
                consider \code{windowSize/2} bases upstream, and \code{windowSize/2-1}
                bases downstream.} For unstranded features, the effect is the same as
                for + strand features.
  \item{type}{What part of the interval to make the window around. If the value is
             \code{"TSS"}, the the start coordinate is used for all + strand features,
             and the end coordinate is used for all - strand features. If \code{"cemter"}
             is chosen, then the coordinate that is half way between the start and end of
             each feature will be used as the reference point. \code{"block"} results in
             the use the start and end coordinates without modification.}
  \item{organism}{The \code{BSgenome} object to calculate mappability upon.}
  \item{verbose}{Whether to print the progess of processing.}
  \item{...}{The \code{verbose} variable for the \code{data.frame} method,
             passed onto the \code{GRanges} method.}
}
\details{
  The windows considered will be \code{windowSize/2} bases upstream and
  \code{windowSize/2-1} bases downstream of the given position of stranded features, 
  and the same number of bases towards the start and end of the chromosome for unstranded
  features. The value returned for each region is a percentage of bases in that region that
  are not N (any base in IUPAC nomenclature).

  For any positions of a window that are off the end of a chromosome, they will be
  considered as being N.
}
\value{
  A vector of mappability percentages, one for each region.
}
\author{Aaron Statham}
\examples{
  \dontrun{
    require(BSgenome.Hsapiens36bp.UCSC.hg18mappability)
    TSSTable <- data.frame(chr = paste("chr", c(1,2), sep = ""), position = c(100000, 200000))
    mappabilityCalc(TSSTable, Hsapiens36bp, 200)
  }
}

