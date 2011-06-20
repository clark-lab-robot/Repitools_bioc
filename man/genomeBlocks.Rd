\name{genomeBlocks}
\alias{genomeBlocks}
\alias{genomeBlocks,numeric-method}
\alias{genomeBlocks,BSgenome-method}
\title{
  Creates bins across a genome.
}
\description{
  Creates a compact \code{GRanges} representation of bins across specified
  chromosomes of a given genome.
}
\usage{
  \S4method{genomeBlocks}{numeric}(genome, chrs = names(genome), width = NULL,
                                   spacing = width)
  \S4method{genomeBlocks}{BSgenome}(genome, chrs = seqnames(genome), width = NULL,
                                    spacing = width)
}
\arguments{
  \item{genome}{Either a \code{BSgenome} object, or a named vector of integers
                (names being choromosome names, integers being the chromosome
                lengths), to get the chromosome lengths from.}
  \item{chrs}{A \code{vector} containing which chromosomes to create bins across.
              May either be numeric indicies or chromosome names. Default is all
              chromosomes given by \code{genome}.}
  \item{width}{The width in base pairs of each bin.}
  \item{spacing}{The space between the centres of each adjacent bin. By default,
                 is equal to the \code{spacing} parameter, which gives
                 non-overlapping bins. Values larger than \code{spacing} will give
                 overlapping bins, and values smaller than  \code{spacing} will
                 give gaps between each bin.}
}
\value{Returns a \code{GRanges} object, compatible with direct usage in
       \code{\link{annotationBlocksCounts}}}
\author{Aaron Statham}
\seealso{\code{\link{annotationBlocksCounts}}}
\examples{
  chr.lengths <- c(800, 200, 200)
  names(chr.lengths) <- c("chr1", "chr2", "chr3")
  genomeBlocks(chr.lengths, width = 200)
}
