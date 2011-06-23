\name{featureScores}
\alias{featureScores}
\alias{featureScores,ANY,data.frame-method}
\alias{featureScores,ANY,GRanges-method}
\title{Get scores at regular sample points around genomic features.}
\description{
  Given a \code{GRanges} / \code{GRangesList} object, or BAM file paths, of reads
  for each experimental condition, or a \code{matrix} or an \code{AffynetrixCelSet},
  and an anntotation of features of interest, scores at regularly spaced positions
  around the features is calculated. In the case of sequencing data, it is the smoothed
  coverage of reads divided by the library size. In the case of array data, it is array
  intensity.
}
\section{Usage}{
  \describe{
    The ANY,data.frame method:
    \code{featureScores(x, anno, ...)}
    The ANY,GRanges method:
    \code{featureScores(x, anno, up = NULL, down = NULL, ...)}
  }
}
\section{Arguments}{
  \describe{
    \item{x:}{If it is a \code{character} vector, then paths to BAM files, otherwise
             a collection of mapped short reads.}
    \item{anno:}{Annotation of the features to sample around.}
    \item{p.anno:}{A \code{data.frame} with columns \code{chr}, \code{position}, an
                  optionally \code{index}. Only provide this if \code{x} is array
                  data. If \code{index} is not provided, the rows are assumed to
                  be in the same order as the elements of \code{x}.}
    \item{mapping:}{A mapping between probes and genes, as made by
                   \code{annotationLookup}. Avoids re-computing the mapping if it
                   has already been done.  Only provide this if \code{x} is array
                   data.}
    \item{up:}{How far to go from the features' reference points in one direction.}
    \item{down:}{How far to go from the features' reference points in the opposite
                direction.}
    \item{dist:}{The type of distance measure to use, in determining the boundaries
                of the sampling area. Only provide this if \code{x} is sequencing
                data. Default: \code{"base"}.}
    \item{freq:}{Score sampling frequency.}
    \item{log2.adj:}{Whether to log2 scale the array intensities. Only provide this
                    if \code{x} is array data. Default: TRUE.}
    \item{s.width:}{The width of smoothing to apply to the coverage. Only provide this
                   if \code{x} is sequencing data. This argument is optional. If
                   not provided, then no smoothing is done.}
    \item{use.strand:}{Whether to only consider reads on the same strand as the
                       feature. Useful for RNA-seq applications.}
    \item{verbose:}{Whether to print the progess of processing. Default: TRUE.}
  }
}
\value{
  A \code{\linkS4class{ScoresList}} object, that holds a list of score
  matrices, one for each experiment type, and the parameters that were used
  to create the score matrices.
}
\details{
  If \code{x} is a vector of paths or \code{GRangesList} object,
  then \code{names(x)} should contain the types of the experiments.

  If \code{anno} is a \code{data.frame}, it must contan the columns \code{chr},
  \code{start}, and \code{end}. Optional columns are \code{strand} and \code{name}.
  If \code{anno} is a \code{GRanges} object, then the name can be present as a column
  called \code{name} in the element metadata of the GRanges object. If names
  are given, then the coverage matrices will use the names as their row names.

  An approximation to running mean smoothing of the coverage is used. Reads are
  extended to the smoothing width, rather than to their fragment size, and
  coverage is used directly. This method is faster than a running mean of the
  calculated coverage, and qualtatively almost identical.

  If providing a matrix of array intensity values, the column names of this
  matrix are used as the names of the samples.

  The annotation can be stranded or not. if the annotation is stranded, then
  the reference point is the start coordinate for features on the + strand,
  and the end coordinate for features on the - strand. If the annotation is
  unstranded (e.g. annotation of CpG islands), then the midpoint of the feature
  is used for the reference point.

  The \code{up} and \code{down} values give how far up and down from the
  reference point to find scores. The semantics of them depend
  on if the annotation is stranded or not. If the annotation is stranded, then
  they give how far upstream and downstream will be sampled. If the annotation is
  unstranded, then \code{up} gives how far towards the start of a chromosome to go,
  and \code{down} gives how far towards the end of a chromosome to go.

  If sequencing data is being analysed, and \code{dist} is \code{"percent"},
  then they give how many percent of each feature's width away from the reference
  point the sampling boundaries are. If \code{dist} is \code{"base"}, then the
  boundaries of the sampling region are a fixed width for every feature, and
  the units of \code{up} and \code{down} are bases. \code{up} and \code{down}
  must be identical if the features are unstranded. The units of \code{freq} are
  percent for \code{dist} being \code{"percent"}, and bases for \code{dist} being
  \code{"base"}.

  In the case of array data, the sequence of positions described by \code{up},
  \code{down}, and \code{freq} actually describe the boundaries of windows, and
  the probe that is closest to the midpoint of each window is chosen as the
  representative score of that window. On the other hand, when analysing sequencing
  data, the sequence of positions refer to the positions that coverage is taken for. 

  Read coverage sampling points that are off the edge of a chromosome (e.g. for a
  gene that has a TSS very near to the edge of a chromosome) are given a coverage
  of 0. Coverage is normalised by dividing the raw coverage by the total number of
  reads in a sample.
}
\seealso{\code{\link{mergeReplicates}} for merging sequencing data replicates of an
         experiment type.}
\author{Dario Strbenac}
\examples{
  data(chr21genes)
  data(samplesList)
  data(expr)

  fs <- featureScores(samplesList[1:2], chr21genes, up = 2000, down = 1000, freq = 500, s.width = 500)
}
