\name{drawChronosCI}
\alias{drawChronosCI}
\title{Draw Confidences Intervals on Phylogenies}
\description{
  This is a low-level plotting command to draw the confidence intervals
  on the node of a tree as rectangles with coloured backgrounds.
}
\usage{
drawChronosCI(CI, col95 = "#FF00004D", col50 = "#0000FF4D", height = NULL,
              legend = TRUE, ...)
}
\arguments{
  \item{CI}{output from \code{\link{chronosCI}} or a similar matrix.}
  \item{col95}{colour used for the 95\% intervals; by default: transparent red.}
  \item{col50}{colour used for the 50\% intervals; by default: transparent blue.}
  \item{height}{the height of the boxes.}
  \item{legend}{a logical value.}
  \item{\dots}{arguments passed to .}
}
\details{
  The matrix \code{CI} must have four rows and as many columns as the
  number of nodes of the tree. The first and fourth rows give the lower
  and upper bounds of the 95\% confidence intervals. The second and
  third rows give the lower and upper bounds of the 50\% confidence
  intervals.
}
\value{NULL}
\author{Emmanuel Paradis, Santiago Claramunt, Joseph Brown, Klaus Schliep}
\seealso{
\code{\link{chronosCI}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
\keyword{aplot}
