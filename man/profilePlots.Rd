\name{profilePlots}
\alias{profilePlots}
\alias{profilePlots,ScoresList-method}
\title{Create line plots of averaged signal across a promoter for gene sets, compared to
       random sampling.}
\description{
  Creates a plot where the average signal across a promoter of supplied gene lists is
  compared to random samplings of all genes, with a shaded confidence area.
}
\usage{
  \S4method{profilePlots}{ScoresList}(x, summarize = c("mean", "median"), gene.lists,
   n.samples = 1000, confidence = 0.975, legend.plot = "topleft", cols = rainbow(length(gene.lists)),
   verbose = TRUE, ...)
}
\arguments{
  \item{x}{A \code{\linkS4class{ScoresList}} object. See \code{\link{featureScores}}.}
  \item{summarize}{How to summarise the scores for each bin into a single value.}
  \item{gene.lists}{Named \code{list} of \code{logical} or \code{integer} vectors, specifying
                    the genes to be averaged and plotted.}
  \item{n.samples}{The number of times to randomly sample from all genes.}
  \item{confidence}{A percentage confidence interval to be plotted (must be > 0.5 and < 1.0).}
  \item{legend.plot}{Where to plot the legend - directly passed to \code{\link{legend}}.
                     NA suppresses the legend.}
  \item{cols}{The colour for each of the genelists supplied.}
  \item{verbose}{Whether to print details of processing.}
  \item{...}{Extra arguments to \code{matplot}, like x- and y-limits, perhaps.}
}
\details{
  For each table of scores in \code{x}, a plot is created showing the average signal of
  the genes specified in each list element of \code{gene.lists} compared to \code{n.samples}
  random samplings of all genes, with \code{confidence} \% intervals shaded. If an element
  of \code{gene.lists} is a \code{logical} vector, its length must be the same as
  the number of rows of the score tables.
}
\value{
  A series of plots.
}
\author{Aaron Statham}
\examples{
  # See examples in manual.
}
