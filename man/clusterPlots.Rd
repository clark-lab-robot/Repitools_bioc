\name{clusterPlots}
\alias{clusterPlots}
\alias{clusterPlots,ClusteredScoresList-method}
\alias{clusterPlots,ScoresList-method}

\title{Visualisation of tables of feature coverages.}
\description{
  Takes the output of \link{featureScores}, or a modified version of it,
  and plots a heatmaps or lineplots representation of clustered coverages.
}
\usage{
  \S4method{clusterPlots}{ClusteredScoresList}(
     c.list, plot.ord = 1:length(c.list), plot.type = c("heatmap", "line", "by.cluster"),
     heat.bg.col = "black", summarize = c("mean", "median"), symm.scale = FALSE, cols = NULL, t.name = NULL,
     verbose = TRUE, ...)
  \S4method{clusterPlots}{ScoresList}(c.list, scale = function(x) x,
    cap.q = 0.95, cap.type = c("sep", "all"), n.clusters = NULL,
    plot.ord = 1:length(c.list), expr = NULL, expr.name = NULL, sort.data = NULL,
    sort.name = NULL, plot.type = c("heatmap", "line", "by.cluster"),
    summarize = c("mean", "median"), cols = NULL, t.name = NULL, verbose = TRUE, ...)
}
\arguments{
  \item{c.list}{A \linkS4class{ScoresList} or \linkS4class{ClusteredScoresList}
                object.}
  \item{scale}{A function to scale all the coverages by. Default : No scaling.}
  \item{cap.q}{The quantile of coverages above which to make any bigger coverages
              equal to the quantile.}
  \item{cap.type}{If \code{"sep"}, then the cap quantile is calculated and applied
                  to each coverage matrix separately. If \code{"all"}, then one cap
                  quantile is calculated based on all of the matrices combined.}
  \item{n.clusters}{Number of clusters to find in the coverage data. Required.}
  \item{plot.ord}{Order of the experiment types to plot.}
  \item{expr}{A vector of expression values.}
  \item{expr.name}{A label, describing the expression data.}
  \item{sort.data}{A vector of values to sort the features within a cluster on.}
  \item{sort.name}{Label to place under the \code{sort.data} plot.}
  \item{plot.type}{Style of plot to draw.}
  \item{heat.bg.col}{If a heatmap is being drawn, the background colour to plot NA
                     values with.}
  \item{summarize}{How to summarise the score columns of each cluster. Not relevant
                  for heatmap plot.}
  \item{symm.scale}{Whether to make lineplot y-axis or heatmap intensity centred
                    around 0. By default, all plots are not symmetrically ranged.}
  \item{cols}{The colours to use for the lines in the lineplot or intensities in
              the heatmap.}
  \item{t.name}{Title to use above all the heatmaps or lineplots. Ignored when 
                cluster-wise lineplots are drawn.}
  \item{verbose}{Whether to print the progress of processing.}
  \item{...}{Further graphical paramters passed to \code{plot} when heatmap plot is
             drawn, that influence how the points of the expression and sort data plots
             will look. If the lineplot is being drawn, parameters to influence
             the line styles.}
}
\value{
  If called with a \code{ScoresList}, then a \code{ClusteredScoresList} is
  returned. If called with a \code{ClusteredScoresList}, then nothing is returned.
}
\details{
  A \code{ClusteredScoresList} should be created by the user, if they wish to do
  some custom clustering and normalisation on the coverage matrices. Otherwise, if
  the user is happy with k-means or PAM clustering, then the \code{ScoresList} object as
  output by \code{featureScores()} can be directly used. If called with a \code{ScoresList},
  then the matrices for each coverage type are joined. Then the function supplied by
  the \code{scale} argument is used to scale the data. Next, each matrix is capped.
  Then each matrix is divided by its maximum value, so that the Euclidean distance
  between maximum reads and no reads is the same for each matrix. Lastly, either k-means
  or PAM clustering is performed to get the cluster membership of each feature. If there are any
  NAs in the scores, then PAM will be used. Otherwise, k-means is used for speed. Then, a
  \code{ClusteredScoresList} object is created, and used. The clusters are
  guaranteed to be given IDs in descending order of summarised cluster expression, if it
  is provided. If called with a \code{ClusteredScoresList}, no scaling or capping
  is done, so it is the user's responsibility to normalise the coverage matrices as
  they see fit, when creating the \code{ClusteredScoresList} object.

  If a \code{ClusteredScoresList} object is subsetted, the original data range is
  saved in a private slot, so that if the user wants to plot a subset of the features,
  such as a certain cluster, for example, the intensity range of the heatmap,
  or the y-axis range of the lineplot will be the same as before subsetting.

  If expression data is given, the summarised expression level of each cluster is
  calculated, and the clusters are plotted in order of decreasing expression, down
  the page. Otherwise, they are plotted in ascending order of cluster ID.
  If a heatmap plot is being drawn, then a heatmap is drawn for every coverage matrix,
  side-by-side, and a plot of each feature's expression is put alongside the heatmaps,
  if provided. If additional sort vector was given, the data within clusters are sorted
  on this vector, then a plot of this data is made as the rightmost graph.

  The lineplot style is similar to the heatmap plot, but clusters are summarised.
  A grid, with as many rows as there are clusters, and as many columns as there are
  clusters is made, and lineplots showing the summarised scores are made in the grid.
  Beside the grid, a boxplot of expression is drawn for each cluster, if provided.

  For a cluster-wise lineplot, a graph is drawn for each cluster, with the colours
  being the different coverage types. Because it makes sense that there will be more
  clusters than there are types of coverage (typically double to triple the number),
  the plots are not drawn side-by-side, as is the layout for the heatmaps. For this
  reason, sending the output to a PDF device is necessary. It is recommended to make
  the width of the PDF device wider than the default. Since the coverage data between
  different marks is not comparable, this method is inappropriate for visualising a
  \code{ClusteredScoresList} object if it was created by the clusterPlots scoresList
  method. If the user, however, can come up with a normalisation method to account
  for the differences that are apparent between different types (i.e. peaked vs.
  spread) of marks that makes the coverages meaningfully comparable, they can alter
  the tables, do their own clustering, and create a \code{ClusteredScoresList}
  object with the modified tables.
}
\seealso{\code{\link{featureScores}} for generating coverage matrices.}
\author{Dario Strbenac}
\examples{
  data(samplesList)
  data(expr)
  data(chr21genes)

  fs <- featureScores(samplesList[1:2], chr21genes, up = 2000, down = 1000, freq = 500, s.width = 500)
  clusterPlots(fs, function(x) sqrt(x), n.clusters = 5, expr = as.numeric(expr),
               plot.type = "heatmap", pch = 19, cex = 0.5)
}
