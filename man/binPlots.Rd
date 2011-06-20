\name{binPlots}
\alias{binPlots}
\alias{binPlots,ScoresList-method}
\title{Create line plots of averaged signal across a promoter}
\description{
  Using a specified ordering of genes, they are split into multiple bins. In each bin,
  the signal across is summarized and displayed visually.
}
\usage{
  \S4method{binPlots}{ScoresList}(x, summarize = c("mean", "median"), ordering,
   ord.label, plot.type = c("line", "heatmap", "terrain"), n.bins = 10, cols = NULL,
   lwd = 3, lty = 1, same.scale = TRUE, symm.scale = FALSE, verbose = TRUE)
}
\arguments{
  \item{x}{A \code{\linkS4class{ScoresList}} object. See \code{\link{featureScores}}.}
  \item{summarize}{How to summarise the scores for each bin into a single value.}
  \item{ordering}{A \code{data.frame} of either numeric or factor variables, with the same
                  number of rows as the annotation used to create \code{x}.}
  \item{ord.label}{Character string that describes what type of data the ordering is.
                   e.g. "log2 expression". Used to label relevant plot axis.}
  \item{plot.type}{Style of plot to draw.}
  \item{n.bins}{The number of bins to split the features into, before summarisation.}
  \item{cols}{A vector of colours to use for the bins. In order from the lowest value
              bin, to the highest value bin.}
  \item{lwd}{Line width of lines in line plot (either scalar or vector).}
  \item{lty}{Line type of line in line plot (either scalar or vector).}
  \item{same.scale}{Whether to keep the scale on all plots be the same.}
  \item{symm.scale}{Whether the scale on plots is symmetrical around 0.}
  \item{verbose}{Whether to print details of processing.}
}
\details{
If \code{plotType = "line"}, a line is plotted for each bin across the promoter.
  
If \code{plotType = "heatmap"}, a series of bins are plotted as a heatmap.  This can be useful to display a larger number of bins.

If \code{plotType = "terrain"}, a series of bins are plotted as a 3D-terrain map.  This can be useful to display a larger number of bins.
}
\value{
  Either a single- or multiple-panel figure.
}
\author{Mark Robinson}
\examples{
  anno.file <- system.file("data", "chr21genes.csv", package = "Repitools")
  annoDF <- read.csv(anno.file)

  reads.file <- system.file("data", "samplesList.RData", package = "Repitools")
  load(reads.file)  # GRangesList of reads 'samplesList'
  
  expr.file <- system.file("data", "expr.RData", package = "Repitools")
  load(expr.file)  # matrix of differential expression 'expr'

  fs <- featureScores(samplesList, annoDF, up = 5000, down = 1000, dist = "base", freq = 1000,
                      s.width = 500)
  fs@scores <- list(tables(fs)[[2]] - tables(fs)[[4]])
  names(fs) <- "PC-Norm"

  binPlots(fs, ordering = expr, ord.label = "expression", plot.type = "line", n.bins = 4)

  binPlots(fs, ordering = expr, ord.label = "expression", plot.type = "heatmap", n.bins = 8)
}
