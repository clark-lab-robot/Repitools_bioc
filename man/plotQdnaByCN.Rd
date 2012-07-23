\name{plotQdnaByCN}
\alias{plotQdnaByCN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting the response of qDNA-seq data by CNV
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotQdnaByCN(obj, cnv.group, idx.ref = 1, idx.sam = 2, min.n = 100, quantile = 0.99, ylim = c(-5, 5), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
%%     ~~Describe \code{obj} here~~
}
  \item{cnv.group}{
%%     ~~Describe \code{cnv.group} here~~
}
  \item{idx.ref}{
%%     ~~Describe \code{idx.ref} here~~
}
  \item{idx.sam}{
%%     ~~Describe \code{idx.sam} here~~
}
  \item{min.n}{
%%     ~~Describe \code{min.n} here~~
}
  \item{quantile}{
%%     ~~Describe \code{quantile} here~~
}
  \item{ylim}{
%%     ~~Describe \code{ylim} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (obj, cnv.group, idx.ref = 1, idx.sam = 2, min.n = 100, 
    quantile = 0.99, ylim = c(-5, 5), ...) 
{
    if (!is(obj, "QdnaData")) 
        stop("Only operates on 'QdnaData' objects.")
    els <- obj$DGEList$samples$lib.size * obj$DGEList$samples$norm.factors
    cnt <- obj$DGEList$counts
    tc <- table(cnv.group)
    n <- length(tc)
    f <- rep(NA, n)
    layout(matrix(c(1:n, rep(n + 1, n)), 2, n, byrow = TRUE))
    for (i in 1:n) {
        cn <- names(tc)[i]
        kk <- cnv.group == cn
        map <- maPlot(cnt[kk, idx.ref]/els[idx.ref], cnt[kk, 
            idx.sam]/els[idx.sam], normalize = FALSE, main = cn, 
            ylim = ylim, ...)
        grid()
        q <- .getQuantile(map, min.n, quantile)
        abline(v = q, col = "blue")
        f[i] <- median(map$M[map$A > q])
        abline(h = f[i], col = "red", lwd = 4)
    }
    f <- exp(f)
    plot(1:n, f, pch = 19, ylab = "Relative Depth by CNV", type = "b", 
        xlab = "", xaxs = "i", xlim = c(1 - 0.5, n + 0.5))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
