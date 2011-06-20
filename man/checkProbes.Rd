\name{checkProbes}
\alias{checkProbes}
\alias{checkProbes,data.frame,data.frame-method}
\alias{checkProbes,GRanges,GRanges-method}
\title{Check Probe Specificity for Some Regions}
\description{Given a set of gene coordinates, and probe mappings to the genome, a plot
  is created across every gene region of how many probes mapped to each position.
}
\usage{
  \S4method{checkProbes}{data.frame,data.frame}(regs, probes, up = NULL, down = NULL, ...)
  \S4method{checkProbes}{GRanges,GRanges}(regs, probes, up = NULL, down = NULL, ...)
}
\arguments{
  \item{regs}{A \code{data.frame} with (at least) columns \code{chr}, \code{start},
              \code{end}, \code{strand}, and \code{name}, or a \code{GRanges} object
              with an elementMetadata column \code{name}. The starts and ends of
              regions describe are the windows plotted in.}
  \item{probes}{A \code{data.frame} describing where the probes mapped to, with
               (at least) columns \code{name} (identifier of a probe), \code{chr},
               \code{start}, and \code{end}, or a \code{GRanges} object with an
               elementMetadata column \code{name}.}
  \item{up}{How many bases upstream to plot.}
  \item{down}{How many bases downstream to plot.}
  \item{...}{Line parameters passed onto \code{matplot}.}
}
\value{
  A set of plots is created, one for each of the genes. The lines
  in the plot show where a probe hits (the x - axis) and how many places in total
  the probe hits in the genome (y - axis).
}
\details{
  If \code{up} and \code{down} are NULL, then the gene is plotted as it is
  described by its start and end coordinates.

  This function produces a number of plots. Sending output to a PDF device is
  recommended.
}
\author{Dario Strbenac}
\examples{
	p.table <- data.frame(name = c("probeA", "probeB", "probeC", "probeC", "probeC"),
			    strand = c('+', '-', '+', '-', '-'),
                               chr = c("chr1", "chr2", "chr1", "chr2", "chr2"),
                             start = c(20, 276, 101, 101, 151),
                               end = c(44, 300, 125, 125, 175))
	r.table <- data.frame(name = c("gene1", "gene2", "gene3"),
                               chr = c("chr1", "chr2", "chr2"),
                            strand = c('+', '-', '+'),
                             start = c(20, 500, 75),
                               end = c(200, 800, 400))
	pdf("tmp.pdf", height = 6, width = 14)
	checkProbes(r.table, p.table, lwd = 4, col = "blue")
	dev.off()
}

