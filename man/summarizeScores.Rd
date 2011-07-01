\name{summarizeScores}
\alias{summarizeScores}
\alias{summarizeScores,ScoresList,matrix-method}
\title{Subtract scores of different samples.}
\description{
  Based on a design matrix, scores matrices are subtracted, and a new \code{ScoresList}
  is returned, with the scores of the contrasts in it.
}
\usage{
  \S4method{summarizeScores}{ScoresList,matrix}(scores.list, design, verbose = TRUE)
}
\arguments{
  \item{scores.list}{A \code{\linkS4class{ScoresList}} object describing the coverage
                     or intensity scores of a set of samples.}
  \item{design}{A matrix that contains only -1, 0, or 1.}
  \item{verbose}{Whether to print a statement explaining the function was called.}
}
\value{
  A \code{\linkS4class{ScoresList}} object holding the scores of the contrasts that were
  specified by the design matrix.
}
\author{Dario Strbenac}
\examples{
  data(chr21genes)
  data(samplesList)
  data(expr)

  fs <- featureScores(samplesList[1:2], chr21genes, up = 2000, down = 1000, freq = 500, s.width = 500)
  d.matrix <- matrix(c(-1, 1))
  colnames(d.matrix) <- "IP-input"
  summarizeScores(fs, d.matrix)
}

