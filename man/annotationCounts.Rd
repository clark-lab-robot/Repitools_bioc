\name{annotationCounts}
\alias{annotationCounts}
\alias{annotationCounts,ANY,data.frame-method}
\alias{annotationCounts,ANY,GRanges-method}
\title{Counts the number of sequencing reads surrounding supplied annotations}
\description{
  Counts are made in windows with boundaries fixed distances either side of a
  reference point.
}
\section{Usage}{
  \describe{
    # ANY,data.frame method
    \code{annotationCounts(x, anno, ...)}
    # ANY,GRanges method
    \code{annotationCounts(x, anno, up, down, ...)}
  }
}
\section{Arguments}{
  \describe{
    \item{x:}{A character vector of BAM paths, \code{GRangesList}, or \code{GRanges} object.}
    \item{anno:}{A set of genomic features to make windows around a reference point of
               theirs. Either a \code{data.frame} with (at least) colums \code{chr},
               \code{start}, and \code{end}, or a \code{GRanges} object.}
    \item{up:}{The number of bases upstream to look.}
    \item{down:}{The number of bases downstream to look.}
    \item{seq.len:}{If sequencing reads need to be extended, the fragment size to be
                   used. Default: NULL (no extension).}
    \item{verbose:}{Whether to print progress. Default: TRUE.}
    \item{...:}{Parameters described above, that are not used in the function called,
               but are passed into \link{annotationBlocksCounts}, that uses them in its
               processing.}
  }
}
\value{A \code{matrix} of counts is returned, one column per sample and one row
       per row of genomic features supplied.}
\details{
    If the genomic features annotation contains all unstranded features, the
    \code{up} and \code{down} distances refer to how far towards the start of
    a chromosome, and how far towards the end to make the counting window
    boundaries. If the annotation is all stranded, then the \code{up} and
    \code{down} distances are relative to the TSS of the features.
}
\author{Aaron Statham}
\seealso{\code{\link{annotationBlocksCounts}}, \code{\link{genomeBlocks}}}
\examples{
  require(GenomicRanges)
  reads <- GRanges(seqnames = rep("chr1", 5),
                   IRanges(c(3309, 4756, 4801, 4804, 5392), width = 36),
                   strand = c('+', '-', '-', '+', '+'))
  genes <- GRanges("chr1", IRanges(5000, 7000), strand = '+')

  annotationCounts(reads, genes, 500, 500, 300)
}
