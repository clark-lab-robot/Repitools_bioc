\name{ScoresList}
\docType{class}
\alias{ScoresList}
\alias{ScoresList-class}
\alias{names,ScoresList-method}
\alias{names<-,ScoresList-method}
\alias{tables}
\alias{tables,ScoresList-method}
\alias{length,ScoresList-method}
\alias{show,ScoresList-method}
\alias{[,ScoresList-method}
\alias{subsetRows}
\alias{subsetRows,ScoresList-method}

\title{Container for \code{featureScores()} output.}
\description{
  Contains a list of tables of sequencing coverages or array intensities,
  and the parameters that were used to generate them.
}
\section{Accessors}{
  In the following code snippets, \code{x} is a ScoresList object.
  
  \describe{
    \item{}{
    	\code{names(x)}, \code{names(x) <- value}
	Gets and sets the experiment type names.
    }
    \item{}{
	\code{tables(x)}
	Gets the list of score matrices.
    }
    \item{}{
	\code{length(x)}
	Gets the number of score matrices.
    }
  }
}
\section{Subsetting}{
  In the following code snippets, \code{x} is a ScoresList object.

  \describe{
    \item{}{
	\code{x[i]}
	Creates a ScoresList object, keeping only the \code{i} matrices.
    }
    \item{}{
	\code{subsetRows(x, i = NULL)}
	Creates a ScoresList object, keeping only the \code{i} features.
    }
  }
}

\author{Dario Strbenac}
