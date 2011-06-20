\name{genQC}
\alias{genQC}
\alias{genQC,character-method}
\alias{genQC,SequenceQCSet-method}

\title{Plot Quality Checking Information for Sequencing Data}
\description{A series of quality control plots for sequencing data are made.}
\usage{
  \S4method{genQC}{character}(qc.data, ...)
  \S4method{genQC}{SequenceQCSet}(qc.data, expt = "Experiment")
}
\arguments{
  \item{qc.data}{A vector of character strings, each containing an absolute path
               to an RData file of a \code{\linkS4class{SequenceQC}} object, or
               a \code{\linkS4class{SequenceQC}} set object.}
  \item{expt}{The names of the experiments which the lanes are about.}
  \item{...}{The \code{expt} argument, which is not directly used in the \code{character}
             method.}
}
\value{
  The function is called for its output. The output is multiple pages, so the
  pdf device should be called before this function is.
}
\details{
  \code{qc.data} can be named, in which case this gives the names of the lanes used
  in the plotting. Otherwise the lanes will be given the names \code{"Lane 1"},
  \code{"Lane 2"}, ..., \code{"Lane n"}.
}
\author{Dario Strbenac}
\references{
  FastQC: \url{http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/}
}
\examples{
  \dontrun{
    qc.files <- list.files(qc.dir, "QC.*RData", full.names = TRUE)
    genQC(qc.files, "My Simple Experiment")
  }
}
