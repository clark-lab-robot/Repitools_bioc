\name{colourBlocks}
\alias{colourBlocks}
\alias{colourBlocks,GRangesList-method}
\title{Plot blocks of features enriched in a particular cluster.}
\description{
  Given a list of blocks, this function plots the anchor position of each feature, with
  width as determined by the window which was used in \code{\link{featureScores}}, on a
  chromosomal illustration.
}
\usage{
  \S4method{colourBlocks}{GRangesList}(blocks, cols = NULL, name.col = NULL, verbose = TRUE)
}
\arguments{
    \item{blocks}{Typically the first element of the list that is output by.}
    \item{cols}{Colours corresponding to the different clusters.}
    \item{name.col}{The column of the metadata of \code{blocks} that contains
                    the feature names.}
    \item{verbose}{Whether to print the progess of processing. Default: TRUE.}
}
\value{
  This function produces a plot for each block. Therefore, output should be sent to PDF.
}

\seealso{\code{\link{closeClusterFeatures}} for how to find blocks.}
\author{Dario Strbenac.}
\examples{
  data(chr21genes)
  data(samplesList)
  data(expr)

  fs <- featureScores(samplesList[1:2], chr21genes, up = 2000, down = 1000, freq = 500, s.width = 500)
  ccl <- clusterPlots(fs, n.clusters = 4, expr = expr, expr.name = "Expression")
  cf <- closeClusterFeatures(ccl, n.consec = 5, n.perm = 10)

  pdf("blocks.pdf", width = 12)
  colourBlocks(cf[[1]], c("red", "green", "blue", "yellow"), 1)
  dev.off()
}



