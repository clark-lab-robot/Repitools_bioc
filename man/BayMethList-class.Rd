\name{BayMethList}
\Rdversion{1.1}
\docType{class}
\alias{BayMethList}
\alias{BayMethList-class}
\alias{BayMethList,GRanges,matrix,matrix,numeric-method}
\alias{[,BayMethList-method}
\alias{length}
\alias{length,BayMethList-method}
\alias{control}
\alias{control<-}
\alias{control,BayMethList-method}
\alias{control<-,BayMethList-method}
\alias{cpgDens}
\alias{cpgDens<-}
\alias{cpgDens,BayMethList-method}
\alias{cpgDens<-,BayMethList-method}
\alias{sampleInterest}
\alias{sampleInterest<-}
\alias{sampleInterest,BayMethList-method}
\alias{sampleInterest<-,BayMethList-method}
\alias{show,BayMethList-method}
\alias{windows}
\alias{windows<-}
\alias{windows,BayMethList-method}
\alias{windows<-,BayMethList-method}
\alias{fOffset}
\alias{fOffset<-}
\alias{fOffset,BayMethList-method}
\alias{fOffset<-,BayMethList-method}
\alias{priorTab}
\alias{priorTab<-}
\alias{priorTab,BayMethList-method}
\alias{priorTab<-,BayMethList-method}
\alias{methEst}
\alias{methEst<-}
\alias{methEst,BayMethList-method}
\alias{methEst<-,BayMethList-method}
\alias{maskEmpBayes}
\alias{maskEmpBayes<-}
\alias{maskEmpBayes,BayMethList-method}
\alias{maskEmpBayes<-,BayMethList-method}
\alias{ncontrol}
\alias{ncontrol,BayMethList-method}
\alias{nsampleInterest}
\alias{nsampleInterest,BayMethList-method}

\title{Class \code{"BayMethList"}}
\description{
    This S4 class captures the genomic windows together with the number of 
    read counts obtained by affinity-enrichment sequencing experiments for a
    fully methylated control and one or more samples of interest. Furthermore
    CpG-density is stored.
}
\section{Constructor}{
    Creates a BayMethList object:
    \describe{
    \item{}{
        \code{BayMethList(windows, control, sampleInterest, cpgDens, 
                f=matrix(), priorTab=list(), methEst=list(), maskEmpBayes=logical())}
        \describe{
            \item{\code{windows}}{A \code{\linkS4class{GRanges}} object.}
            \item{\code{control}}{A matrix of read counts obtained by an
                affinity enrichment sequencing experiment for the fully
                methylated (SssI) treated sample. The number of rows must be
                equal to \code{length(windows)}. Each column contains the
                counts of one sample. The number of columns must be either one
                or equal to the number of columns of \code{sampleInterest}.}
            \item{\code{sampleInterest}}{A matrix of read counts obtained by an
                affinity enrichment sequencing experiment for the samples of
                interest. The number of rows must be equal to
                \code{length(windows)}. Each column contains the counts of
                one sample. }
            \item{\code{cpgDens}}{A numeric vector containing the CpG density
                for \code{windows}. The length must be equal to
                \code{length(windows)}}
            \item{\code{fOffset}}{A matrix where each column contains
                the normalizing offsets for one sample. The number of rows must
                be either equal to one or the number of windows.}
            \item{\code{priorTab}}{A list containing for each sample of interest
                the prior parameters as determined by \code{empBayes}.}
            \item{\code{methEst}}{A list containing the methylation estimates
                as determined by \code{methylEst}.}
            \item{\code{maskEmpBayes}}{A logical vector indicating which bins should
                be masked out in the empirical Bayes analysis. TRUE indicates to 
                neglect the bin in the empirical Bayes approach.}
        }
    }
  }
}
\section{Methods}{
  \describe{
    \item{x[i]}{\code{signature(x = "BayMethList")}: Creates a BayMethList 
        object, keeping only the \code{i} entries. }
    \item{length}{\code{signature(x= "BayMethList")}: gets the number of 
        genomic regions included. }
    \item{control<-}{\code{signature(x = "BayMethList")}: replace the 
        \code{control} slot}
    \item{control}{\code{signature(object = "BayMethList")}: extract the 
        \code{control} matrix slot. }
    \item{cpgDens<-}{\code{signature(x = "BayMethList")}: replace the
        \code{cpgDens} slot}
    \item{cpgDens}{\code{signature(object = "BayMethList")}:  extract the
        \code{cpgDens} slot. }
    \item{sampleInterest<-}{\code{signature(x = "BayMethList")}: 
        replace the \code{sampleInterest} slot}
    \item{sampleInterest}{\code{signature(object = "BayMethList")}: 
        extract the \code{sampleInterest} matrix slot.}
    \item{show}{\code{signature(object = "BayMethList")}: 
        show an overview of the object }
    \item{windows<-}{\code{signature(x = "BayMethList")}: 
        replace the \code{windows} slot }
    \item{windows}{\code{signature(object = "BayMethList")}: 
        extract the \code{windows} GRanges slot. }
    \item{fOffset<-}{\code{signature(x = "BayMethList")}: 
        replace the \code{fOffset} slot }
    \item{fOffset}{\code{signature(object = "BayMethList")}: 
        extract the \code{fOffset} slot. }
    \item{priorTab<-}{\code{signature(x = "BayMethList")}: 
        replace the \code{priorTab} slot }
    \item{priorTab}{\code{signature(object = "BayMethList")}: 
        extract the \code{priorTab} slot. }
    \item{methEst<-}{\code{signature(x = "BayMethList")}: 
        replace the \code{methEst} slot }
    \item{methEst}{\code{signature(object = "BayMethList")}: 
        extract the \code{methEst} slot. }
    \item{maskEmpBayes<-}{\code{signature(x = "BayMethList")}: 
        replace the \code{maskEmpBayes} slot }
    \item{maskEmpBayes}{\code{signature(object = "BayMethList")}: 
        extract the \code{maskEmpBayes} slot. }
    \item{ncontrol}{\code{signature(object = "BayMethList")}: 
        get the number of provided SssI samples.}
    \item{nsampleInterest}{\code{signature(object = "BayMethList")}: 
        get the number of provided samples of Interest.}
	 }
}

\seealso{determineOffset, empBayes, methylEst }

\author{
    Andrea Riebler and Mark Robinson
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

        cat("Number of genomic regions", length(bm), "\n")
        cat("Number of fully methylated control samples:", ncontrol(bm), "\n")
        cat("Number of samples of interest:", nsampleInterest(bm), "\n")
        bm[2:20]
    }
}
\keyword{classes}
