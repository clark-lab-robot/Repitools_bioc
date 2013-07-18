\name{determineOffset}
\alias{determineOffset}
\title{
    Function to determine the normalising offset f that accounts 
    for the relative sequencing depth.
}
\description{
The composition of a library influences the resulting read densities. 
To adjust the modelled mean (in the Poisson model) for these composition
effects, we estimate a normalising factor f that accounts simultaneously for
overall sequencing depth and composition. The derivation of this offset is
motivated by the M (log ratio) versus A (average-log-count) plot. 
}
\usage{
determineOffset(x, quantile = 0.998, controlPlot = list(show = FALSE, 
    nsamp = 50000, mfrow=c(1,1), xlim=NULL, ylim=NULL, main=NULL, ask=FALSE))
}

\arguments{
  \item{x}{
    \code{BayMethList} object.
}
  \item{quantile}{
    quantile q to restrict values of A = log2(sampleInterest*control)/2 
}
  \item{controlPlot}{
    list defining whether a MA plot should be shown. 
    \describe{
    \item{-}{
        \code{show}
        logical. If 'TRUE' the corresponding MA plot is shown. (default FALSE)
    }
    \item{-}{
        \code{nsamp}
        number of genomic regions included in the plot. 
        (These are sampled without replacement).
    }
    \item{-}{
        \code{mfrow}
        vector of the form "c(nr, nc)" to determine how several plots 
        should be ordered.
    }
    \item{-}{
        \code{xlim, ylim}
        numeric vectors of length 2, giving the x and y coordinates ranges.
    }
    \item{-}{
        \code{main}
        If NULL the names of the sample of interest are used as title in the MA plot. 
        Alternatively, a vector with length equal to the number of samples of interest
        can be provided.
    }
    \item{-}{
        \code{ask}
        logical.  If 'TRUE' (and the R session is interactive) the
          user is asked for input, before a new figure is drawn. (default FALSE).
        }
    }   
  }
}
\value{
    A \code{BayMethList} object given as input, where the slot \code{fOffset}
    is filled accordingly.
}
\author{
    Andrea Riebler
}
\seealso{
    maPlot, plotSmear
}
\examples{
   if(require(BSgenome.Hsapiens.UCSC.hg18)){
        windows <- genomeBlocks(Hsapiens, chrs="chr21", width=100, spacing=100)
        cpgdens <- cpgDensityCalc(windows, organism=Hsapiens, 
            w.function="linear", window=700)  
        co <- matrix(rnbinom(length(windows), mu=10, size=2), ncol=1)
        sI <- matrix(rnbinom(2*length(windows), mu=5, size=2), ncol=2)
        bm <- BayMethList(windows=windows, control=co, sampleInterest=sI,
            cpgDens=cpgdens)

        bm <- determineOffset(bm, controlPlot=list(show=TRUE, mfrow=c(1,2)))
   }
}

