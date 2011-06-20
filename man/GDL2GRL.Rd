\name{GDL2GRL}
\alias{GDL2GRL}
\alias{GDL2GRL,GenomeDataList-method}
\title{Utility function to covert a GenomeDataList object into GRangesList objects.}
\description{The data in the \code{\linkS4class{GenomeDataList}} object is made into
             a \code{\linkS4class{GRangesList}} object.}
\usage{
  \S4method{GDL2GRL}{GenomeDataList}(gdl)
}
\arguments{
  \item{gdl}{A GenomeDataList.}
}
\value{
  A \code{\linkS4class{GRangesList}}.
}
\author{Dario Strbenac}
\examples{
    require(BSgenome)
    gdl <- GenomeDataList(list(
                         GenomeData(list(
                                         chr1 = list(`-` = c(100, 200), `+` = c(800, 1000)),
                                         chr2 = list(`-` = c(450, 550), `+` = c(1500, 7500))
                                         )
                                    ),
                         GenomeData(list(
                                         chr1 = list(`-` = c(300, 700), `+` = c(850, 900)),
                                         chr2 = list(`-` = c(125, 250), `+` = c(500, 750))
                                         )
                                    )
                                )
                            )
    GDL2GRL(gdl)
}
