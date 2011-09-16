\name{featureScores}
\alias{featureScores}
\alias{featureScores,ANY,data.frame-method}
\alias{featureScores,ANY,GRanges-method}
\title{Get scores at regular sample points around genomic features.}
\description{
  Given a \code{GRanges} / \code{GRangesList} object, or BAM file paths, of reads
  for each experimental condition, or a \code{matrix} or an \code{AffynetrixCelSet},
  or a numeric matrix of array data, where the rows are probes and the columns are
  the different samples,and an anntotation of features of interest, scores at
  regularly spaced positions around the features is calculated. In the case of
  sequencing data, it is the smoothed coverage of reads divided by the library size.
  In the case of array data, it is array intensity.
}
\section{Usage}{
  \describe{
    The ANY,data.frame method: \cr
    \code{featureScores(x, anno, ...)} \cr
    The ANY,GRanges method: \cr
    \code{featureScores(x, anno, up = NULL, down = NULL, ...)}
  }
}
\section{Arguments}{
  \describe{
    \item{x:}{Paths to BAM files, a collection of mapped short reads, or a collection of
              microarray data.}
    \item{anno:}{Annotation of the features to sample around.}
    \item{p.anno:}{A \code{data.frame} with columns \code{chr}, \code{position}, an
                  optionally \code{index}. Only provide this if \code{x} is array
                  data. If \code{index} is not provided, the rows are assumed to
                  be in the same order as the elements of \code{x}.}
    \item{mapping:}{A mapping between probes and genes, as made by
                   \code{annotationLookup}. Avoids re-computing the mapping if it
                   has already been done.  Only provide this if \code{x} is array
                   data.}
    \item{chrs:}{A mapping between chromosome names in an ACP file to the user's
                 feature annotation. Only provide this if \code{x} is an
                 \code{AffymetrixCelSet}. There is no need to provide this if
                 the feature annotation uses the same chromosome names as the ACP
                 files do. Element \code{i} of this vector is the name to give to
                 the chromosome numbered \code{i} in the ACP information.}
    \item{up:}{How far to go from the features' reference points in one direction.}
    \item{down:}{How far to go from the features' reference points in the opposite
                direction.}
    \item{dist:}{The type of distance measure to use, in determining the boundaries
                of the sampling area. Only provide this if \code{x} is sequencing
                        data. Default: \code{"base"}. \code{"percent"} is also accepted.}
    \item{freq:}{Score sampling frequency.}
    \item{log2.adj:}{Whether to log2 scale the array intensities. Only provide this
                    if \code{x} is array data. Default: TRUE.}
    \item{s.width:}{The width of smoothing to apply to the coverage. Only provide this
                   if \code{x} is sequencing data. This argument is optional. If
                   not provided, then no smoothing is done.}
    \item{mappability:}{A \code{\linkS4class{BSgenome}} object that has bases for which
                        no mappable reads start at masked by N. If this was provided, then
                        either \code{s.width} or \code{tag.len} must be provided (but not
                        both).}
    \item{map.cutoff:}{The percentage of bases in a window around each sampling position
                       that must be mappable. Otherwise, the score at that position is
                       repalced by NA. Default: 0.5}
    \item{tag.len:}{Provide this if \code{mappability} was provided, but \code{s.width}
                    was not.}
    \item{use.strand:}{Whether to only consider reads on the same strand as the
                       feature. Useful for RNA-seq applications.}
    \item{verbose:}{Whether to print the progess of processing. Default: TRUE.}
  }
}

\section{Details}{
  \describe{
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

  Providing a mappability object for sequencing data is recommended. Otherwise, it is
  not possible to know if a score of 0 is because the window around the sampling position
  is unmappable, or if there were really no reads mapping there in the experiment.
  Coverage is normalised by dividing the raw coverage by the total number of
  reads in a sample. The coverage at a sampling position is multiplied by 1 / mappability.
  Any positions that have mappabilty below the mappability cutoff will have their
  score set to NA.
}}

\section{Value}{
  \describe{
    A \code{\linkS4class{ScoresList}} object, that holds a list of score
    matrices, one for each experiment type, and the parameters that were used
    to create the score matrices.
  }
}

\seealso{\code{\link{mergeReplicates}} for merging sequencing data replicates of an
         experiment type.}
\author{Dario Strbenac, with contributions from Matthew Young at WEHI.}
\examples{
  data(chr21genes)
  data(samplesList)
  data(expr)

  fs <- featureScores(samplesList[1:2], chr21genes, up = 2000, down = 1000, freq = 500, s.width = 500)
}
