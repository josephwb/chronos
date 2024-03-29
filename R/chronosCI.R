## chronosCI.R (2021-10-28)
## Copyright 2021 Emmanuel Paradis, Santiago Claramunt, Joseph Brown, Klaus Schliep

##' @title Confidence Intervals of Chronograms
##' @description Computes the confidence intervals of a chronogram using
##' different methods.
##' @param chronogram a chronogram output by [ape::chronos()].
##' @param pml.output an unrooted tree output by [phangorn::optim.pml()].
##' @param B the number of replications.
##' @param type the bootstrap method to be used. Can be either a character
##' string among "nonparametric", "semiparametric", or "parametric", or an
##' integer between 1 and 3.
##' @param calibration a data frame with the calibration points of the tree.
##' @param method a character string specifying the method used to assess
##' calibration uncertainty. See [qage()] for available methods.
##' @param control see [ape::chronos()]. By default, these are the same than
##' used to estimate the chronogram with [chronos()].
##' @param quiet a logical value. By default, the progress of the
##' computations is printed.
##' @param trees a logical value specifying whether to return the bootstrap
##' trees produced by phangorn (FALSE, by default).
##' @details The details of the methods are presented in the manuscript below.
##'
##' The labels (or taxa names) of the first argument (\code{chronogram})
##' must be present in the second argument (\code{pml.output}), and this
##' second one must also include the outgroup.
##' @return a matrix with the 50\% and 95\% lower and upper bounds of the
##' confidence intervals.
##' @references Paradis, E., Claramunt, S., Brown, J., and Schliep, K. Confidence
##' intervals in molecular dating by penalized likelihood. (in preparation)
##' @author Emmanuel Paradis, Santiago Claramunt, Joseph Brown, Klaus Schliep
##' @importFrom ape chronos.control Nedge.phylo drop.tip root makeNodeLabel Ntip branching.times chronos
##' @importFrom stats rpois
##' @importFrom phangorn bootstrap.pml pml.control
##' @importFrom stats quantile
##' @examples
##' \dontrun{
##' ##---- Should be DIRECTLY executable !! ----
##' ##-- ==>  Define data, use random,
##' ##--	or do  help(data=index)for the standard data sets.
##' }
##' @seealso [qage()] [drawChronosCI()] [ape::chronos()]
##' [phangorn::bootstrap.pml()]
##' @keywords models
##' @export
chronosCI <-
    function(chronogram, pml.output, B = 100, type = "semiparametric",
             calibration = NULL, method = "StraussSadler",
             control = NULL, quiet = FALSE, trees = FALSE)
{
    getArgument <- function(name, default) {
        i <- match(name, args)
        if (is.na(i)) return(default)
        res <- thecall[[i]]
        if (!is.character(res)) eval(res) else res
    }
    thecall <- attr(chronogram, "call")
    args <- names(thecall) # the arguments used in the call
    model <- getArgument("model", "correlated")
    lambda <- getArgument("lambda", 1)
    if (is.null(control)) control <- getArgument("control", chronos.control())

    BOOTSTRAPMETHODS <- c("nonparametric", "semiparametric", "parametric")
    if (is.character(type))
        type <- pmatch(type, BOOTSTRAPMETHODS)
    if (is.na(type) || type > length(BOOTSTRAPMETHODS) || !is.numeric(type))
        stop("bootstrap method cannot be found")
    switch(type, {
        smoothRates <- FALSE
        Poisson <- FALSE
    }, {
        smoothRates <- TRUE
        Poisson <- FALSE
    }, {
        smoothRates <- FALSE
        Poisson <- TRUE
    })

    treeML <- pml.output$tree
    taxa <- treeML$tip.label
    ingroup <- chronogram$tip.label
    outgroup <- taxa[! taxa %in% ingroup]
    if (length(outgroup) < 1 || all(is.na(outgroup)))
        stop("no outgroup identified from the arguments")

    if (Poisson) {
        if (!quiet) cat("Generating Poisson branch lengths...")
        S <- length(attr(pml.output$data, "index"))
        tr0 <- treeML
        tr0$tip.label <- NULL
        TR <- vector("list", B)
        N <- ape::Nedge.phylo(treeML)
        for (i in 1:B) {
            tr0$edge.length <- stats::rpois(N, S * treeML$edge.length) / S
            TR[[i]] <- tr0
        }
        class(TR) <- "multiPhylo"
        attr(TR, "TipLabel") <- taxa
    } else {
        if (!quiet) cat("Running phangorn bootstrap...")
        TR <- phangorn::bootstrap.pml(pml.output, B, optNni = FALSE,
                            control = phangorn::pml.control(trace = 0L))
    }
    if (!quiet) cat("\n")

    ## root and drop the outgroup from the bootstrap trees:
    rTR <- lapply(TR, function(x) ape::drop.tip(ape::root(x, outgroup), outgroup))

    ## make labels associated to the edges of the chronogram:
    tmp <- c(ingroup, ape::makeNodeLabel(chronogram, "md5sum")$node.label)
    original.edge.labels <- tmp[chronogram$edge[, 2]]
    original.ints <- chronogram$edge[, 2] > ape::Ntip(chronogram)
    ## id. for the 1st bootstrap tree:
    tr1 <- rTR[[1]]
    tmp <- c(tr1$tip.label, ape::makeNodeLabel(tr1, "md5sum")$node.label)
    boot.edge.labels <- tmp[tr1$edge[, 2]]
    boot.ints <- tr1$edge[, 2] > ape::Ntip(tr1)
    ## to reorder the results below
    needReorder <- !identical(original.edge.labels, boot.edge.labels)
    if (needReorder)
        o <- match(original.edge.labels[original.ints], boot.edge.labels[boot.ints])
    ## a finir....

    ## get bootstrap branch lengths:
    if (smoothRates && !Poisson) {
        EL <- sapply(rTR, "[[", "edge.length")
        NE <- nrow(EL)
        bootBL <- matrix(NA_real_, NE, B)
        for (i in 1:NE) {
            f <- getDensel(EL[i, ])
            F <- Cdensel(f)
            bootBL[i, ] <- rbl(B, F)
        }
        for (i in 1:B) rTR[[i]]$edge.length <- bootBL[, i]
    }

    ##
    CALS <- list()
    if (!is.null(calibration)) {
        calNodes <- unique(calibration$node)
        for (nd in calNodes) {
            i <- which(calibration$node == nd)
            if (length(i) == 1) next
            ages <- as.matrix(calibration[i, c("age.min", "age.max")])
            cals <- rage(ages, n = B, method = method)
            CALS <- c(CALS, list(cals))
            names(CALS)[length(CALS)] <- nd
        }
    } else {
        cal <- data.frame(node = length(rTR[[1]]$tip.label) + 1L,
                          age.min = 1, age.max = 1, soft.bounds = FALSE)
    }

    CHR <- vector("list", B)

    if (length(CALS)) {
        nodes2 <- as.integer(names(CALS))
        DUP <- duplicated.default(calibration$node)
    } else {
        cal <- calibration
    }

    for (i in 1:B) {
        if (!quiet) cat("\rRunning chronos bootstrap:", i, "/", B)
        if (length(CALS)) {
            cal <- calibration[!DUP, ]
            ages2 <- sapply(CALS, "[", i)
            cal[match(names(ages2), cal$node), 2:3] <- ages2
        }
        CHR[[i]] <- ape::chronos(rTR[[i]], quiet = TRUE, model = model,
                                 calibration = cal, control = control)
    }
    BT <- sapply(CHR, ape::branching.times)
    CI <- apply(BT, 1, stats::quantile, probs = c(0.025, 0.25, 0.75, 0.975))
    if (!quiet) cat("\nDone.\n")
    if (needReorder) CI <- CI[, o]
    if (!trees) return(CI)
    class(rTR) <- "multiPhylo"
    list(CI = CI, trees = rTR)
}

isUnimodal <- function(density)
{
    v <- rle(sign(diff(density$y)))$values
    if (length(v) > 3) return(FALSE)
    if (length(v) == 3 && all(v == c(1, 0, -1))) return(TRUE)
    if (length(v) == 2 && all(v == c(1, -1))) return(TRUE)
    FALSE
}

getDensel <- function(el)
{
    start <- stats::bw.nrd0(el)
    inc <- start / 10
    b <- start
    repeat {
        densel <- stats::density(el, bw = b)
        if (isUnimodal(densel)) {
            xsmall <- densel$x < 1e-8
            if (any(xsmall)) {
                w <- which(xsmall)
                densel$x <- densel$x[-w]
                total <- sum(densel$y)
                extra <- sum(densel$y[w])
                densel$y <- densel$y[-w]
                densel$y <- densel$y * total / (total - extra)
            }
            break
        }
        b <- b + inc
    }
    densel
}

Cdensel <- function(f)
{
    n <- length(f$y)
    list(x = (f$x[-1] + f$x[-n]) / 2,
         y = cumsum(diff(f$x) * (f$y[-1] + f$y[-n]) / 2))
}

## random branch lengths
rbl <- function(n, F) {
    z <- numeric(n)
    p <- stats::runif(n)
    N <- length(F$y)
    if (any(s <- p < F$y[1])) z[s] <- F$x[1]
    if (any(l <- p > F$y[N])) z[l] <- F$x[N]
    m <- which(!s & !l)
    if (!length(m)) return(z)
    for (k in m) {
        j <- which(F$y > p[k])[1]
        i <- j - 1L
        b <- (F$y[j] - F$y[i]) / (F$x[j] - F$x[i])
        a <- F$y[i] - b * F$x[i]
        z[k] <- (p[k] - a) / b
    }
    z
}
