\name{annotationBlocksLookup}
\alias{annotationBlocksLookup}
\alias{annotationBlocksLookup,data.frame,data.frame-method}
\alias{annotationBlocksLookup,data.frame,GRanges-method}
\title{Forms a mapping between probe locations and chromosomal blocks (regions).}
\description{
  Starting from a table of genome locations for probes, and a table of regions
  of interest, this procedure forms a list structure that contains the indices
  to map from one to the other.
}
\usage{
  \S4method{annotationBlocksLookup}{data.frame,data.frame}(x, anno, ...)
  \S4method{annotationBlocksLookup}{data.frame,GRanges}(x, anno, verbose = TRUE)
}
\arguments{
  \item{x}{probe genomic locations, a \code{data.frame} with required elements
           \code{chr}, \code{position}, and optionally \code{index}}
  \item{anno}{a \code{data.frame} with required elements \code{chr}, \code{start},
              \code{end}, \code{strand} and optional element \code{name}.
              Also may be a \code{GRanges} with optional elementMetadata column
              \code{name}.}
  \item{verbose}{Whether to print progress to screen.}
  \item{...}{Represents the \code{verbose} parameter, when the \code{data.frame,data.frame}
             method is called.}
}
\details{
  Strandedness of probes is ignored, even if it is given.

  If \code{x} has no index column, then the probes are given indices from 1 to
  the number of probes, in the order that they appear in the \code{data.frame} or
  \code{GRanges} object. 
}
\value{
  A list with elements
  \item{indexes}{a list for each gene in \code{y}, giving a vector
                 of indices to the probe data.}
  \item{offsets}{a list for each gebe in \code{y}, giving a vector
                 (corresponding to \code{indexes}) of offsets relative to the
                 start of the block.}
}
\author{Aaron Statham, Mark Robinson}
\seealso{\code{\link{annotationLookup}} which simplifies annotation lookups for
        constant sized regions}
\examples{
# create example set of probes and gene start sites
probeTab <- data.frame(position=seq(1000,3000,by=200), chr="chrX", strand="+")
genes <- data.frame(chr="chrX", start=c(2100,2200), end=c(2500, 2400), strand=c("+","-"))
rownames(genes) <- paste("gene",1:2,sep="")

# Call annotationLookup() and look at output
annotationBlocksLookup(probeTab, genes)
}
