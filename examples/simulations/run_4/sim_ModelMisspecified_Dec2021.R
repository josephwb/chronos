## Objectives: assess confidence intervals when the model is wrong

source("../COMMON_CODE.R")
library(chronos)

nrep <- 100L # number of replications
## rate.sim <- 1 # rate used by phangorn::simSeq()
pmlctr <- pml.control(trace = 0L)

## best to not erase outputs from the previous simulations:
OUTFILE <- getOutfile("sim_wrongmodel", "rda")
## NOTE: don't put the end-of-prefix number (001, 002, ...) in the above call

## the parameters to vary
## NTIPS <- c(20, 50, 100) # number of species (taxa, sequences)
NCALPT <- c(1, 5, 10)
## S <- c(1e3, 1e4) # number of nucleotides in the sequences
MODEL <- c(1, 2) # the studied models

mod4chronos <- c("clock", "correlated")

## create a data frame with the parameter combinations:
PARA <- expand.grid(MODEL = MODEL, NCALPT = NCALPT)

## NOTE: this way, all simulations are run with two nested loops: the
## first one running along the rows of PARA, and the second one
## replicating the simulations for each parameter combinations.

OUT <- vector("list", nrep * nrow(PARA))

set.seed(3) # seed used for the simulations reported in the paper
## modify or delete the previous line to have different results

n <- 20L
s <- 1000L

k <- 0L
for (i in 1:nrow(PARA)) {
    Ncalpt <- PARA$NCALPT[i]
    model <- PARA$MODEL[i]

    ## np <- n + Ntip(outtree)

    for (j in 1:nrep) {
        k <- k + 1L
        cat(k, file = "k")
        tr <- simTreeWithOutgroup(n)
        trsim <- rescaleEdgeLength(FUNTREESIM[[model]](tr))

        cal <- getCalibrationPoints(Ncalpt)
        X <- simData(trsim, s, 1)

        ## ML estimation with phangorn
        utr <- unroot(trsim)
        op.out <- optim.pml(pml(utr, X, rate = rate.sim), control = pmlctr)
        tr.ml <- root(op.out$tree, "outgroup")
        tr4chr <- drop.tip(tr.ml, "outgroup")

        if (Ncalpt > 1) {
            for (i in 2:Ncalpt) {
                nodetr <- cal$node[i]
                NODE <- getMRCA(tr4chr, PP[[nodetr - n]])
                if (NODE != nodetr) cal$node[i] <- NODE
            }
        }

        wrongmodel <- rev(mod4chronos)[model]

        lambda <- 1
        if (wrongmodel == "correlated")
            lambda <-  getLambdaCV(tr4chr, wrongmodel, cal)

        chr <- chronos(tr4chr, lambda, wrongmodel, !TRUE, cal)

        ## nonparametric boostrap
        res1 <- try(chronosCI(chr, op.out, 100, 1, cal, quiet=!TRUE), silent=TRUE)
        if (class(res1) == "try-error") res1 <- matrix(Inf, 4, n - 1)
        ## semiparametric boostrap
        res2 <- try(chronosCI(chr, op.out, 100, 2, cal, quiet=!TRUE), silent=TRUE)
        if (class(res2) == "try-error") res2 <- matrix(Inf, 4, n - 1)
        ## parametric boostrap
        res3 <- try(chronosCI(chr, op.out, 100, 3, cal, quiet=!TRUE), silent=TRUE)
        if (class(res3) == "try-error") res3 <- matrix(Inf, 4, n - 1)
        BT <- branching.times(drop.tip(tr, "outgroup"))
        OUT[[k]] <- cbind(BT, t(res1), t(res2), t(res3))
        ## save results at each loop
        save(OUT, file = OUTFILE)
    }
}
