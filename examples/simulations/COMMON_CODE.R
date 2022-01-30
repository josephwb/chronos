## NOTE: NELSI has some dependencies not defined in its DESCRIPTION file;
## to avoid a full installation of NELSI with all dependencies (not all available on CRAN),
## unzip the source (NELSI-master.zip) and souce() the files:
setwd("NELSI_RCodeOnly/")
source("get.tree.data.matrix.R")
source("mid.edge.ages.R")
source("simulate.clock.R")
source("simulate.uncor.lnorm.R")
source("simulate.uncor.gamma.R")
source("simulate.autocor.kishino.R")
source("allnode.times.R")
setwd("../")
## library(NELSI)
library(phangorn)

FUNTREESIM <-
list(function(phy) simulate.clock(phy, params = list(rate = 0.01, noise = 1e-5))[[1]],
function(phy) simulate.autocor.kishino(phy, params = list(initial.rate = 1e-5, v = 0.3))[[1]],
function(phy) simulate.uncor.lnorm(phy, params = list(mean.log = log(0.01), sd.log = 0.1))[[1]],
function(phy) simulate.uncor.gamma(phy, list(mean = NULL, shape = 0.5, rate = 1))[[1]])
##function(phy) simulate.uncor.gamma(phy, list(mean = 0.001, shape = 3.98, rate = NULL))[[1]])

##simBaseTree <- function(n, threshold = 1, limit = 1e4) {
##    sperate <- switch(as.character(n),
##                      "20" = 0.1, "50" = 0.05, "100" = 0.025, "200" = 0.01)
##    for (l in seq_len(limit)) {
##        tr <- rphylo(n, sperate, 0)
##        ## rescale so that the root age is 50:
##        tr$edge.length <- 50 * tr$edge.length/branching.times(tr)[1]
##        return(tr)
##        terms <- tr$edge[, 2] <= n
##        ##terms[] <- TRUE
##        if (all(tr$edge.length[terms] > threshold)) return(tr)
##    }
##    warning("reached the limit of ", limit, " tries: returned the tree anyway")
##    tr
##}

simBaseTree <- function(n) {
    SPERATE <- c("20" = 0.1, "50" = 0.05, "100" = 0.025, "200" = 0.01)
    tr <- rphylo(n, SPERATE[as.character(n)], 0)
    ## rescale so that the root age is 50:
    tr$edge.length <- 50 * tr$edge.length / branching.times(tr)[1]
    tr
}

outtree <- read.tree(text = "(outgroup:1);")
simTreeWithOutgroup <- function(n, root.edge = 1) {
    tr <- simBaseTree(n)
    ## look at branching times before grafting the outgroup
    BTtr <<- branching.times(tr)
    PP <<- prop.part(tr) # to get MRCAs later
    mrca.age <<- BTtr[1]
    ## add the outgroup
    tr$root.edge <- root.edge
    otree <- outtree
    otree$edge.length <- mrca.age + root.edge
    tr + otree
}

## create data drame with calibration times
## (mrca.age and BTtr are found by lexical scoping)
getCalibrationPoints <- function(Ncalpt = 1) {
    cal <- data.frame(node = n + 1L, age.min = mrca.age,
                      age.max = mrca.age, soft.bounds = FALSE)
    if (Ncalpt > 1) {
        extra <- sample.int(n - 2L, Ncalpt - 1L) + 1L
        ## extra varies between 2 and n - 1
        tmp <- BTtr[extra]
        DF <- data.frame(node = extra + n, age.min = tmp,
                         age.max = tmp, soft.bounds = FALSE)
        cal <- rbind(cal, DF)
    }
    cal
}

## rescale branch lengths so that the shortest terminal branch length
## is at least threshold (0.001 by default)
rescaleEdgeLength <- function(phy, threshold = 0.001) {
    shortest.terminal <- min(phy$edge.length[phy$edge[, 2] > Ntip(phy)])
    if (shortest.terminal < threshold)
        phy$edge.length <- phy$edge.length * threshold / shortest.terminal
    phy
}
## rational behind the rescaling:
## Q <- matrix(1e-3/3, 4, 4); diag(Q) <- 0; diag(Q) <- -rowSums(Q)
## matexpo(Q)[1]^s

## threshold: maximum proportion of polymorphic sites
simData <- function(phy, s, rate = 1, threshold = 0.5) {
    repeat {
        X <- simSeq(phy, l = s, rate = rate)
        if (attr(X, "nr") <= threshold * s) break
        rate <- rate / 1.25
    }
    rate.sim <<- rate
    X
}

getOutfile <- function(prefix, fileext = "out", keep.old = TRUE,
                       digits = 3, silent = FALSE)
{
    OUTFILE <- paste(prefix, fileext, sep = ".")
    if (file.exists(OUTFILE)) {
        if (keep.old) {
            regexp <- paste0("^", prefix, ".*\\.", fileext, "$")
            fls <- grep(regexp, dir(), value = TRUE)
            Ntmp <- 1L
            regexp2 <- paste0("[[:digit:]]{", digits, "}\\.", fileext)
            if (length(flsd <- grep(regexp2, fls, value = TRUE))) {
                from <- nchar(prefix) + 1L
                to <- from + digits - 1L
                Ntmp <- 1L + max(as.numeric(substr(flsd, from, to)))
            }
            Ntmp <- sprintf(paste0("%0", digits, "d"), Ntmp)
            newOUTFILE <- paste0(prefix, Ntmp, ".", fileext)
            if (!silent) message(paste("\n\tA file named", sQuote(OUTFILE), "is in the working directory;\n\tusing", sQuote(newOUTFILE), "instead.\n"))
            OUTFILE <- newOUTFILE
        } else {
            if (!silent) message(paste("\n\tA file", sQuote(OUTFILE), "already exists;\n\tit will be overwritten.\n"))
        }
    }
    OUTFILE
}

getLambdaCV <- function(phy, model, calibration) {
    cv <- CV(tr4chr, 10^(-3:3), model, TRUE, calibration)
    cv[which.min(cv[, 2]), 1]
}
