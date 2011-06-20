\name{blocksStats}
\alias{blocksStats}
\alias{blocksStats,ANY,data.frame-method}
\alias{blocksStats,ANY,GRanges-method}

\title{Calculate statistics for regions in the genome}
\description{
  For each region of interest or TSS, this routine interrogates probes or sequence
  data for either a high level of absolute signal or a change in signal for some
  specified contrast of interest.  Regions can be surroundings of TSSs, or can be
  user-specified regions. The function determines if the \code{start} and \code{end}
  coordinates of \code{anno} should be used as regions or as TSSs, if the up and down
  coordinates are \code{NULL} or are numbers.
}
\section{Usage}{
  \describe{
    The ANY,data.frame method:
    \code{blocksStats{ANY,data.frame}(x, anno, ...)}
    The ANY,GRanges method:
    \code{blocksStats{ANY,GRanges}(x, anno, up = NULL, down = NULL, ...)}
  }
}
\section{Arguments}{
  \describe{
    \item{x:}{A \code{GRangesList}, \code{AffymetrixCelSet},
             or a \code{data.frame} of data. Or a \code{character} vector of BAM
             paths to the location of the BAM files.}
    \item{anno:}{Either a \code{data.frame} or a \code{GRanges} giving the gene
                coordinates or regions of interest. If it is a \code{data.frame},
                then the column names are (at least) \code{chr}, \code{name},
                \code{start}, \code{end}. Column \code{strand} is also mandatory,
                if \code{up} and \code{down} are \code{NULL}.}
    \item{seq.len:}{If sequencing reads need to be extended, the fragment size to be used.}
    \item{p.anno:}{A \code{data.frame} with (at least) columns \code{chr},
                  \code{position}, and \code{index}. This is an optional parameter of
                  the \code{AffymetrixCelSet} method, because it can be automatically
                  retrieved for such array data. The parameter is also optional, if
                  \code{mapping} is not \code{NULL}.}
    \item{mapping:}{If a mapping with \code{annotationLookup} or \code{annotationBlocksLookup}
                   has already been done, it can be passed in, and avoids unnecessary
                   re-conmputing of the mapping list within \code{blocksStats}.}
    \item{chrs:}{If \code{p.anno} is \code{NULL}, and is retrieved from an ACP file, this
                 vector gives the textual names of the chromosomes.}
    \item{log2.adj:}{Whether to take $log_2$ of array intensities.}
    \item{design:}{A design matrix specifying the contrast to compute (i.e. The samples
                  to use and what differences to take.).}
    \item{up:}{The number of bases upstream to consider in calculation of statistics.
              If not provided, the starts and ends in \code{anno} are used as
              region boundaries.}
    \item{down:}{The number of bases upstream to consider in calculation of statistics.
                If not provided, the starts and ends in \code{anno} are used as
                region boundaries.}
    \item{lib.size:}{A string that indicates whether to use the total lane count, 
                    total count within regions specified by \code{anno}, or
                    normalisation to a reference lane by the negative binomial
                    quantile-to-quantile method, as the library size for each lane.
                    For total lane count use \code{"lane"}, for region sums use
                    \code{"blocks"}, and for the normalisation use \code{"ref"}.}
    \item{robust:}{Numeric. If it is 0, then a robust linear model is not fitted.
                  If it is greater than 0, a robust linear model is used, and the
                  number specifies the minimum number of probes a region has to
                  have, for statistics to be reported for that region.}
    \item{p.adj:}{The method used to adjust p-values for multiple testing. Possible
                 values are listed in \code{\link{p.adjust}}.}
    \item{Acutoff:}{If \code{libSize} is \code{"ref"}, this argument must be provided.
                   Otherwise, it must not. This parameter is a cutoff on the "A"
                   values to take, before calculating trimmed mean.}
    \item{verbose:}{Logical; whether to output commments of the processing.}
    \item{...}{Parameters described above, that are not used in the function called,
               but are passed further into a private function that uses them in its
               processing.}
  }
}
\details{
  For array data, the statstics are either determined by a t-test, or a linear
  model. For sequencing data, the two groups are assumed to be from a
  negative binomial distribution, and an exact test is used.
}
\value{
  A \code{data.frame}, with the same number of rows as there are features described
  by \code{anno}, but with additional columns for the statistics calculated at each
  feature.
}
\author{Mark Robinson}
\seealso{\code{\link{annotationLookup}} and \code{\link{annotationBlocksLookup}}}
\examples{
  require(GenomicRanges)
  intensities <- matrix(c(6.8, 6.5, 6.7, 6.7, 6.9,
                          8.8, 9.0, 9.1, 8.0, 8.9), ncol = 2)
  colnames(intensities) <- c("Normal", "Cancer")
  d.matrix <- matrix(c(-1, 1))
  colnames(d.matrix) <- "Cancer-Normal"
  probe.anno <- data.frame(chr = rep("chr1", 5),
                           position = c(4000, 5100, 6000, 7000, 8000), 
                           index = 1:5)
  anno <- GRanges("chr1", IRanges(7500, 10000), '+', name = "Gene 1")
  blocksStats(intensities, anno, 2500, 2500, probe.anno, log2.adj = FALSE, design = d.matrix)
}
