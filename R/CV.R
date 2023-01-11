## simple version of drop.tip(tr, tip, collapse.singles = FALSE)
.droptip <- function(phy, tip) {
    s <- which(phy$edge[, 2] == tip)
    phy$edge.length <- phy$edge.length[-s]
    phy$edge <- phy$edge[-s, ]
    k <- which(phy$edge > tip)
    phy$edge[k] <- phy$edge[k] - 1L
    phy$tip.label <- phy$tip.label[-tip]
    phy
}

##' @title Cross-Validation with Penalized Likelihood Chronograms
##' @description Computes the cross-validation (CV) of a chronogram
##'     using penalized likelihood.
##' @param phy an object of class \code{"phylo"}.
##' @param lambda a numeric vector with the values of smoothing
##'     parameter to be assessed.
##' @param model the model of rate evolution. currently set to autocorrelated.
##' @param quiet a logical value.
##' @param calibration a data frame with the calibration dates of the
##'     chronogram.
##' @details This follows the CV calculation described in Sanderson (2002).
##' @return a matrix with two columns named \code{lambda} and
##'     \code{CV}.
##' @references Paradis, E., Claramunt, S., Brown, J., and Schliep,
##'     K. Confidence intervals in molecular dating by penalized
##'     likelihood. (in preparation)
##' @references Sanderson, M. J. (2002) Estimating absolute rates of
##'     molecular evolution and divergence times: a penalized
##'     likelihood approach. \emph{Molecular Biology and Evolution},
##'     \bold{19}, 101--109.
##' @author Emmanuel Paradis, Santiago Claramunt, Joseph Brown, Klaus Schliep
##' @importFrom ape Ntip branching.times chronos makeChronosCalib
##' @examples
##' ## a simple random tree, so lambda is expected to be zero
##' library(ape)
##' tr <- rtree(10)
##' res <- CV(tr, 10^seq(-4, 2, 0.25))
##' plot(res, type = "o", log = "xy")
##'
##' ## create 'artificial' auto-correlation among branches by sorting
##' ## them according to their lengths:
##' tr$edge.length <- sort(tr$edge.length)
##' res <- CV(tr, 10^seq(-4, 2, 0.25))
##' plot(res, type = "o", log = "xy")
##' @seealso [qage()] [drawChronosCI()] [ape::chronos()] [ape::makeChronosCalib()]
##' @keywords models
##' @export
CV <- function(phy, lambda = 10^(-3:2), model = "correlated",
               quiet = FALSE, calibration = makeChronosCalib(phy))
{
    cv <- numeric()
    n <- Ntip(phy)
    nb <- n - 1L
    D <- numeric(n)
    one2n <- 1:n
    S <- match(one2n, phy$edge[, 2])
    cal2 <- calibration
    cal2$node <- cal2$node - 1L
    for (l in lambda) {
        if (!quiet)
            cat("Running cross-validation with lambda = ", l, "\n", sep = "")
        chr <- chronos(phy, l, model, quiet = TRUE, calibration = calibration)
        for (i in one2n) {
            if (!quiet)
                cat("\r    Dropping terminal branch ", i, "/", n, sep = "")
            trb <- .droptip(phy, i)
            chrb <- chronos(trb, l, model, quiet = TRUE, calibration = cal2)
            s <- S[i]
            ancb <- phy$edge[s, 1L] - 1L
            ## ancb is the node# of the direct ancestor of the dropped tip *in* the modified tree
            btb <- branching.times(chrb)
            j <- if (ancb == nb + 1L) 1L else 2L
            xstar <- btb[ancb - nb] * attr(chrb, "rates")[which(trb$edge[, j] == ancb)]
            D[i] <- (xstar - phy$edge.length[s])^2 / xstar
        }
        sD <- sum(D)
        if (!quiet) cat("   ==>   CV = ", sD, "\n", sep = "")
        cv <- c(cv, sD)
    }
    cbind(lambda = lambda, CV = cv)
}
