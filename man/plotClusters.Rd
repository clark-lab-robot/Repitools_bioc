\name{plotClusters}
\alias{plotClusters}
\alias{plotClusters,data.frame-method}
\alias{plotClusters,GRanges-method}
\title{Plot Scores of Cluster Regions}
\description{Given an annotation of gene positions that has a score column, the
             function will make a series of bar chart plots, one for each cluster.
}
\usage{
  \S4method{plotClusters}{data.frame}(x, s.col = NULL, non.cl = NULL, ...)
  \S4method{plotClusters}{GRanges}(x, s.col = NULL, non.cl = NULL, ...)
}
\arguments{
  \item{x}{A summary of genes and their statistical score, and the cluster that
           they belong to. Either a \code{data.frame} or a \code{GRanges}.
           If a \code{data.frame}, then (at least) columns \code{chr}, \code{start},
           \code{end}, \code{strand}, \code{name} and \code{cluster}. Also a score
           column, with the column name describing what type of score it is. If a
           \code{GRanges}, then the \code{elementMetadata} should have a
           \code{DataFrame} with a score column, and columns named \code{"cluster"}
           and \code{"name"}.}
  \item{s.col}{The column number of the \code{data.frame} when \code{data} is a
               \code{data.frame}, or the column number of the \code{DataFrame}
               when \code{data} is a \code{GRanges} object. The name of this
               column is used as the y-axis label in the plot.}
  \item{non.cl}{The value in the cluster column that represents genes that
                are not in any cluster}
  \item{...}{Further parameters to be passed onto \code{plot}.}
}
\value{
  A plot for each cluster is made. Therefore, the PDF device should be opened
  before this function is called.
}
\author{Dario Strbenac}
\examples{
  library(GenomicRanges)
  g.summary <- GRanges("chr1",
                       IRanges(seq(1000, 10000, 1000), width = 100),
                       rep(c('+', '-'), 5),
                       `t-statistic` = rnorm(10, 8, 2),
                       cluster = c(0, 0, 0, 0, 0, 1, 1, 1, 1, 0),
                       name = paste("Gene", 1:10))
  plotClusters(g.summary, 1, 0, ylim = c(4, 12), lwd = 5)
}
