\name{featureBlocks}
\alias{featureBlocks}
\alias{featureBlocks,data.frame-method}
\alias{featureBlocks,GRanges-method}
\title{Make windows for distances around a reference point.}
\description{
  Windows are made around a reference point, which is the start coordinate for
  features on the + strand, and the end coordinate for features on the - strand.
  For unstranded features, the reference point is taken to be the mid-point of
  the feature.
}
\usage{
  \S4method{featureBlocks}{data.frame}(anno, ...)
  \S4method{featureBlocks}{GRanges}(anno, up = NULL, down = NULL, dist = c("base", "percent"),
                                    keep.strand = FALSE)
}
\arguments{
  \item{anno}{A \code{data.frame} or \code{GRanges}, describing some genomic features.}
  \item{up}{The amount to go upstream or towards the start of a chromosome.
            Semantics depend on the value of \code{dist}. See details.}
  \item{down}{The amount to go downstream or towards the end of a chromosome.
            Semantics depend on the value of \code{dist}. See details.}
  \item{dist}{Whether \code{up} and \code{down} refer to bases, or a percentage
              of each feature's width.}
  \item{keep.strand}{Whether the blocks should keep the strands of their features,
                     or if all blocks should have strand be \code{'*'}}
  \item{...}{Arguments from the list above that are not used directly within the
             \code{data.frame} method.}
}
\details{
  \code{up} refers to how many bases to go upstream for stranded features, or
  for unstranded features, how many bases to go towards the start of the chromosome,
  from the mid-point of the feature. Having a negative value for \code{up} means
  that the windows will start downstream by that amount, for stranded features.
  For unstranded features, it will start that many bases closer to the end of
  the chromosome, relative to the feature mid-point.

  \code{down} is defined analogously.
}
\value{
  A \code{\linkS4class{GRanges}} of windows surrounding reference points for the
  features described by \code{anno}.
}
\author{Dario Strbenac}
\examples{
  genes <- data.frame(chr = c("chr1", "chr3", "chr7", "chr22"),
                   start = seq(1000, 4000, 1000),
                   end = seq(1500, 4500, 1000),
                   strand = c('+', '-', '-', '+'))

  featureBlocks(genes, 500, 500)
}

