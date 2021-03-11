## chronosCI.R (2021-03-11)
## Copyright 2021 Emmanuel Paradis, Santiago Claramunt, Joseph Brown, Klaus Schliep

## This file is part of the R-package `chronos'.
## See the file ../LICENSE for licensing issues.

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
    start <- bw.nrd0(el)
    inc <- start / 10
    b <- start
    repeat {
        densel <- density(el, bw = b)
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
    p <- runif(n)
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

chronosCI <-
    function(chronogram, pml.output, B = 100, smoothRates = TRUE,
             Poisson = FALSE, calibration = NULL, method = "StraussSadler",
             quiet = FALSE)
{
    call.str <- deparse(attr(chronogram, "call"))
    if (length(grep("model", call.str))) {
        model <- gsub("^.*model {0,}= {0,}\"", "", call.str)
        model <- gsub("\".*$", "", model)
    } else {
        model <- "correlated"
    }

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
        N <- Nedge.phylo(treeML)
        for (i in 1:B) {
            tr0$edge.length <- rpois(N, S * treeML$edge.length) / S
            TR[[i]] <- tr0
        }
        class(TR) <- "multiPhylo"
        attr(TR, "TipLabel") <- taxa
    } else {
        if (!quiet) cat("Running phangorn bootstrap...")
        TR <- bootstrap.pml(pml.output, B, optNni = FALSE,
                            control = pml.control(trace = 0L))
    }
    if (!quiet) cat("\n")

    ## root and drop the outgroup from the bootstrap trees:
    rTR <- lapply(TR, function(x) drop.tip(root(x, outgroup), outgroup))

    ## make labels associated to the edges of the chronogram:
    tmp <- c(ingroup, makeNodeLabel(chronogram, "md5sum")$node.label)
    original.edge.labels <- tmp[chronogram$edge[, 2]]
    original.ints <- chronogram$edge[, 2] > Ntip(chronogram)
    ## id. for the 1st bootstrap tree:
    tr1 <- rTR[[1]]
    tmp <- c(tr1$tip.label, makeNodeLabel(tr1, "md5sum")$node.label)
    boot.edge.labels <- tmp[tr1$edge[, 2]]
    boot.ints <- tr1$edge[, 2] > Ntip(tr1)
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
        CHR[[i]] <- chronos(rTR[[i]], quiet = TRUE, model = model, calibration = cal)
    }
    BT <- sapply(CHR, branching.times)
    CI <- apply(BT, 1, quantile, probs = c(0.025, 0.25, 0.75, 0.975))
    if (!quiet) cat("\nDone.\n")
    if (needReorder) CI <- CI[, o]
    CI
}

drawChronosCI <- function(CI, col95 = "#FF00004D", col50 = "#0000FF4D", height = NULL,
                          legend = TRUE, ...)
{
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    present <- lastPP$xx[1]
    if (is.null(height)) height <- yinch(0.2)

    ## the nodes are renumbered => need to use labels

    for (i in 1:ncol(CI)) {
        L <- present - CI[1, i]
        R <- present - CI[4, i]
        B <- lastPP$yy[i + lastPP$Ntip] - height / 2
        T <- B + height
        rect(L, B, R, T, col = col95, border = NULL)
        L <- present - CI[2, i]
        R <- present - CI[3, i]
        rect(L, B, R, T, col = col50, border = NULL)
    }
    if (!identical(legend, FALSE)) {
        loc <- if (is.logical(legend)) "topleft" else legend
        legend(loc, legend = c("95% CI", "50% CI"), pch = 22,
               pt.bg = c(col95, col50), col = c(col95, col50), ...)
    }
}

