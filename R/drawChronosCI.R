
##' @title Draw Confidences Intervals on Phylogenies
##' @description This is a low-level plotting command to draw the confidence
##' intervals on the node of a tree as rectangles with coloured backgrounds.
##' @param CI output from [chronosCI()] or a similar matrix.
##' @param col95 colour used for the 95% intervals; by default: transparent
##' red.
##' @param col50 colour used for the 50% intervals; by default: transparent
##' blue.
##' @param height the height of the boxes.
##' @param legend a logical value.
##' @param \dots arguments passed to [graphics::legend()]
##' @details The matrix \code{CI} must have four rows and as many columns as the
##' number of nodes of the tree. The first and fourth rows give the lower and
##' upper bounds of the 95% confidence intervals. The second and third rows
##' give the lower and upper bounds of the 50% confidence intervals.
##' @return NULL
##' @author Emmanuel Paradis, Santiago Claramunt, Joseph Brown, Klaus Schliep
##' @importFrom graphics legend rect yinch
##' @examples 
##' \dontrun{
##' ##---- Should be DIRECTLY executable !! ----
##' ##-- ==>  Define data, use random,
##' ##--	or do  help(data=index)for the standard data sets.
##' }
##' @seealso [chronosCI()]
##' @keywords aplot
##' @export
drawChronosCI <- function(CI, col95 = "#FF00004D", col50 = "#0000FF4D",
                          height = NULL, legend = TRUE, ...)
{
  lastPP <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  present <- lastPP$xx[1]
  if (is.null(height)) height <- graphics::yinch(0.2)
  
  ## the nodes are renumbered => need to use labels
  
  for (i in 1:ncol(CI)) {
    L <- present - CI[1, i]
    R <- present - CI[4, i]
    B <- lastPP$yy[i + lastPP$Ntip] - height / 2
    T <- B + height
    graphics::rect(L, B, R, T, col = col95, border = NULL)
    L <- present - CI[2, i]
    R <- present - CI[3, i]
    graphics::rect(L, B, R, T, col = col50, border = NULL)
  }
  if (!identical(legend, FALSE)) {
    loc <- if (is.logical(legend)) "topleft" else legend
    graphics::legend(loc, legend = c("95% CI", "50% CI"), pch = 22,
           pt.bg = c(col95, col50), col = c(col95, col50), ...)
  }
}
