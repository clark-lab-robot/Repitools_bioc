\name{ClusteredScoresList}
\docType{class}
\alias{ClusteredScoresList}
\alias{ClusteredScoresList-class}
\alias{ClusteredScoresList,ScoresList-method}
\alias{show,ClusteredScoresList-method}
\alias{[,ClusteredScoresList-method}
\alias{subsetRows,ClusteredScoresList-method}
\alias{clusters}
\alias{clusters,ClusteredScoresList-method}

\title{Container for coverage matrices with clustering results.}
\description{
  Contains a list of coverage matrices, the parameters that were used to
  generate them origin, and also cluster membership and expression data.

  This allows \code{\link{clusterPlots}} to be called multiple times without
  having to redo the clustering each time. It also allows the user to take
  the \code{\linkS4class{ScoresList}} output of \code{\link{featureScores}},
  and do their own custom clustering on the coverage matrices, then save the
  clustering results in this container.
}
\section{Constructor}{
  \describe{
    \item{}{
    	\code{ClusteredScoresList(x, cvgs = tables(x), expr = NULL, expr.name = NULL,
                                    cluster.id, sort.name = NULL, sort.data = NULL)}
	Creates a ClusteredScoresList object.
	\describe{
	    \item{\code{x}}{A \code{\linkS4class{ScoresList}} object.}
	    \item{\code{scores}}{A list of coverage matrices.}
	    \item{\code{expr}}{A numeric vector, same length as number of
                               rows of every coverage matrix.}
	    \item{\code{expr.name}}{A label, describing the expression data.}
	    \item{\code{cluster.id}}{A vector, same length as number of
                                     rows of every coverage matrix.}
	    \item{\code{sort.data}}{Vector of data to order features
                                     within clusters by.}
	    \item{\code{sort.name}}{Human readable description of what
                                     the sorting data is of.}
	}
    }
  }
}
\section{Subsetting}{
  In the following code snippets, \code{x} is a ClusteredScoresList object.

  \describe{
    \item{}{
	\code{x[i]}
	Creates a ClusteredScoresList object, keeping only the \code{i} matrices.
    }
    \item{}{
	\code{subsetRows(x, i = NULL)}
	Creates a ClusteredScoresList object, keeping only the \code{i} features.
    }
    \item{}{
	\code{clusters(x)}
	Creates a ClusteredScoresList object, keeping only the \code{i} features.
    }
  }
}
\section{Accessors}{
  In the following code snippets, \code{x} is a ClusteredScoresList object.

  \describe{
    \item{}{
	\code{clusters(x)}
	Get the cluster ID of each feature.
    }
  }
}

\author{Dario Strbenac}
