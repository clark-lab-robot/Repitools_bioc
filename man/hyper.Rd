\name{hyperg2F1_vec}
\alias{hyperg2F1_vec}
\title{Gaussian hypergeometric function for vectorial arguments}
\description{
    Computes the value of the Gaussian hypergeometric function 2_F_1 as 
    defined in Abramowitz and Stegun (1972, page 558), i.e. for |z| < 1 and 
    c > b > 0 using the Cephes library.
}
\usage{
hyperg2F1_vec(a,b,c,z)
}
\arguments{
  \item{a}{(Vectorial) parameter a.}
  \item{b}{parameter b (of same length as a)}
  \item{c}{parameter c (of same length as a)}
  \item{z}{parameter z (of same length as a)}
}
\details{
The function is in particular efficient for vectorial arguments as the 
loop is shifted to \code{C}. Note: If vectorial arguments are provided, all
arguments need to be of the same length.
}
\value{
The value of the Gaussian hypergeometric function F(a,b,c,z) for c > b > 0 and
|z| < 1.
}
\references{
Abramowitz and Stegun 1972.  _Handbook of mathematical functions
with formulas, graphs and mathematical tables_.  New York: Dowver
Publications.

www.netlib.org/cephes/
}
\author{Andrea Riebler and Daniel Sabanes Bove}
\seealso{package \code{hypergeo} or \code{BMS}.}
\examples{
    hyperg2F1_vec(-10.34,2.05,3.05,0.1725)
    hyperg2F1_vec(30,1,20,.8) # returns about 165.8197
    hyperg2F1_vec(30,10,20,0) # returns one
    hyperg2F1_vec(10,15,20,-0.1) # returns about 0.4872972

    hyperg2F1_vec(c(-10.34, 30, 10), c(2.05, 1, 10), c(3.05, 20, 20), 
        c(0.1725, 0.8, 0))
    hyperg2F1_vec(a=1.2+1:10/10, b=rep(1.4,10), c=rep(1.665,10), z=rep(.3,10))
}
\keyword{ programming }
