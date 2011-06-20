\name{annotationLookup}
\alias{annotationLookup}
\alias{annotationLookup,data.frame,data.frame-method}
\alias{annotationLookup,data.frame,GRanges-method}
\title{Forms a mapping between probes on a tiling array and windows surrounding the
       TSSs of genes.}
\description{
  Starting from genome locations for probes and a locations for a set of genes,
  this procedure forms a list structure that contains the indices to map from one
  to the other.
}
\section{Usage}{
  \describe{
    The data.frame,data.frame method:
    \code{annotationLookup(x, anno, ...)}
    The data.frame,GRanges method:
    \code{annotationLookup(x, anno, up, down, ...)}
  }
}
\section{Arguments}{
  \describe{
    \item{x:}{Probe genomic locations, a \code{data.frame} with required elements
           \code{chr}, \code{position}, and optionally \code{index}}
    \item{anno:}{a \code{data.frame} with required elements \code{chr}, \code{start},
               \code{end}, \code{strand} and optional element \code{name}.
               Also may be a \code{GRanges} with optional elementMetadata column
               \code{name}.}
    \item{up:}{The number of bases upstream to look.}
    \item{down:}{The number of bases downstream to look.}
    \item{verbose:}{Whether to print progress to screen. Default: TRUE}
    \item{...:}{Parameters described above, that are not used in the function called,
               but are passed further into \link{annotationBlocksLookup}, which uses
               them in its processing.}
  }
}
\details{
  This function is a wrapper for the generic function \code{annotationBlocksLookup}
  which can handle annotations of varying sizes. \code{annotationLookup} is
  appropriate where you wish to map probes that are within a fixed distance of
  points of annotation e.g gene transcription start sites. Even if strand information
  is given for probes, it is ignored.

  If \code{x} has no index column, then the probes are given indices from 1 to
  the number of probes, in the order that they appear in the \code{data.frame} or
  \code{GRanges} object.

  It is an error for the gene annotation to have unstranded features.
}
\value{
  A list with elements
  \item{indexes}{a list for each gene in \code{y}, giving a vector
                 of indices to the probe data.}
  \item{offsets}{a list for each gebe in \code{y}, giving a vector
                 (corresponding to \code{indexes}) of offsets relative to the
                 genes' TSSs for each probe that mapped that that gene.}
}
\author{Aaron Statham, Mark Robinson}
\seealso{\code{\link{annotationBlocksLookup}}, \code{\link{makeWindowLookupTable}}}
\examples{

# create example set of probes and gene start sites
probes <- data.frame(position=seq(1000, 3000, by = 200), chr = "chrX", strand = '-')
genes <- data.frame(chr = "chrX", start=c(2100, 1000), end = c(3000, 2200),
                    strand=c("+","-"))
rownames(genes) <- paste("gene", 1:2, sep = '')

# Call annotationLookup() and look at output
annotationLookup(probes, genes, 500, 500)
}
