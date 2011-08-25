\name{closeClusterFeatures}
\alias{closeClusterFeatures}
\alias{closeClusterFeatures,ClusteredScoresList-method}
\title{Find features in clusters that are co-located.}
\description{
  Given a \code{\linkS4class{ClusteredScoresList}} object, the function finds
  genomic regions where the features from a cluster occur close together. A
  randomisation procedure is done, and a p-value calculated to give
  the statistical significance of the blocks of features.
}
\usage{
  \S4method{closeClusterFeatures}{ClusteredScoresList}(c.list, n.consec = NULL,
                                 n.perm = 100, name.col = NULL, verbose = TRUE)
}
\arguments{
    \item{c.list}{A \code{\linkS4class{ClusteredScoresList}} object.}
    \item{n.consec}{Each block reported must have at least this many features in a
                     consecutive sequence from the one cluster. This argument is
                      mandatory.}
    \item{n.perm}{The number of times to re-order the cluster IDs and perform block
                   finding.}
    \item{name.col}{The column of the element metadata of \code{c.list} which has the
                    feature IDs. This argument is optional.}
    \item{verbose}{Whether to print the progess of processing. Default: TRUE.}
}

\details{
  The algorithm works by first finding all genomic regions where a cluster is present
  in disperse way. These disperse blocks have features of a cluster that are not
  separated by more than one feature from a different cluster. The disperse blocksS are
  only kept if there is a run of features belonging to the cluster at least \code{n.consec}
  long.

  The p-value is calculated as the number of randomisations that have the same or
  more blocks found than the experimental data does.
}

\value{
    A \code{\linkS4class{ClusteredScoresList}} object, that holds a list of score
    matrices, one for each experiment type, and the parameters that were used
    to create the score matrices.
}

\seealso{\code{\linkS4class{ClusteredScoresList}} for how to create the appropriate object.}
\author{Dario Strbenac.}
\examples{
  data(chr21genes)
  data(samplesList)
  data(expr)

  fs <- featureScores(samplesList[1:2], chr21genes, up = 2000, down = 1000, freq = 500, s.width = 500)
  ccl <- clusterPlots(fs, n.clusters = 4, expr = expr, expr.name = "Expression")
  closeClusterFeatures(ccl, n.consec = 5, n.perm = 10, name.col = 2)
}

