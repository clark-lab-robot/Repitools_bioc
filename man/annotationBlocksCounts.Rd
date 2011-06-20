\name{annotationBlocksCounts}
\alias{annotationBlocksCounts}
\alias{annotationBlocksCounts,ANY,data.frame-method}
\alias{annotationBlocksCounts,character,GRanges-method}
\alias{annotationBlocksCounts,GRanges,GRanges-method}
\alias{annotationBlocksCounts,GRangesList,GRanges-method}
\title{Counts the number of sequencing reads within supplied genomic blocks.}
\description{
  Counts reads inside blocks.
}
\usage{
  \S4method{annotationBlocksCounts}{ANY,data.frame}(x, anno, ...)
  \S4method{annotationBlocksCounts}{character,GRanges}(x, anno, ...)
  \S4method{annotationBlocksCounts}{GRanges,GRanges}(x, anno, seq.len = NULL, verbose = TRUE)
  \S4method{annotationBlocksCounts}{GRangesList,GRanges}(x, anno, ...)
}
\arguments{
  \item{x}{A character vector of BAM paths, a \code{GRangesList}, or \code{GRanges} object.}
  \item{anno}{A set of genomic features to make windows around a reference point of
             theirs. Either a \code{data.frame} with (at least) colums \code{chr},
             \code{start}, and \code{end}, or a \code{GRanges} object.}
  \item{seq.len}{If sequencing reads need to be extended, the fragment size to be
                 used. Default: NULL (no extension).}
  \item{verbose}{Whether to print progress. Default: TRUE.}
  \item{...}{Parameters described above, that are not used in the top-level
             error-checking stage, but are passed further into a private
             function that uses them in its processing.}
}
\value{
  A \code{matrix} of counts is returned, one column per sample and one row
  per row of genomic features supplied.
}
\author{Aaron Statham}
\seealso{\code{\link{annotationCounts}}, \code{\link{genomeBlocks}}}
\examples{
  require(GenomicRanges)
  reads <- GRanges(seqnames = rep("chr1", 5),
                   IRanges(c(3309, 4756, 4801, 4804, 5392), width = 36),
                   strand = c('+', '-', '-', '+', '+'))
  genes <- GRanges("chr1", IRanges(5000, 7000), strand = '+')
  annotationBlocksCounts(reads, genes, 300)
}
