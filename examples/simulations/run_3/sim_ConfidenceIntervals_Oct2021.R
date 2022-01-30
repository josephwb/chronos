## Objectives: assess confidence intervals

source("../COMMON_CODE.R")
library(chronos)
## library(pegas) # to check unique sequences below

nrep <- 100L # number of replications
##rate.sim <- 1 # rate used by phangorn::simSeq(); set later by simData()
# chrctr <- chronos.control(dual.iter.max = 100)
pmlctr <- pml.control(trace = 0L)

## best to not erase outputs from the previous simulations:
OUTFILE <- getOutfile("sim_CI", "rda")

## the parameters to vary (or not if vector of length one)
NTIPS <- 20 #c(20, 50, 100) # number of species (taxa, sequences)
NCALPT <- c(1, 10) #c(1, 5, 10) # number of calibration points
S <- c(1e3, 1e4) # number of nucleotides in the sequences
MODEL <- c(1, 2, 4) # the models

mod4chronos <- c("clock", "correlated", "relaxed", "relaxed")

## create a data frame with the parameter combinations:
PARA <- expand.grid(MODEL = MODEL, NCALPT = NCALPT, S = S, NTIPS = NTIPS)

## NOTE: all simulations are run with two nested loops: the first one
## running along the rows of PARA, and the second one replicating the
## simulations for each parameter combinations.

OUT <- vector("list", nrep * nrow(PARA))
RUNNINGTIMES <- matrix(0, nrep * nrow(PARA), 3)

set.seed(2) # seed used for the simulations reported in the paper
## modify or delete the previous line to have different results

k <- 0L
for (i in 1:nrow(PARA)) {
    n <- PARA$NTIPS[i]
    Ncalpt <- PARA$NCALPT[i]
    s <- PARA$S[i]
    model <- PARA$MODEL[i]

##    np <- n + Ntip(outtree) # needed if pegas::haplotype is called

    for (j in 1:nrep) {
        k <- k + 1L
        cat(k, file = "k")
        tr <- simTreeWithOutgroup(n)
        trsim <- rescaleEdgeLength(FUNTREESIM[[model]](tr))
## graphically compare both trees:
## layout(matrix(1:2, 2))
## plot(tr); axisPhylo()
## plot(trsim); axisPhylo()

        cal <- getCalibrationPoints(Ncalpt)
        X <- simData(trsim, s, 1)
## check number of uniques sequences:
## hap <- haplotype(as.DNAbin(X))
## if (nrow(hap) < np) warning("some sequences were identical")

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

        lambda <- 1
        if (model != 1)
            lambda <-  getLambdaCV(tr4chr, mod4chronos[model], cal)

        chr <- chronos(tr4chr, lambda, mod4chronos[model], TRUE, cal)#, control = chrctr)

        t0 <- proc.time()[3]
        ## nonparametric boostrap
        res1 <- try(chronosCI(chr, op.out, 100, 1, cal, quiet=TRUE), silent=TRUE)
        t1 <- proc.time()[3]
        RUNNINGTIMES[k, 1] <- t1 - t0
        if (class(res1) == "try-error") res1 <- matrix(Inf, 4, n - 1)
        t0 <- proc.time()[3]
        ## semiparametric boostrap
        res2 <- try(chronosCI(chr, op.out, 100, 2, cal, quiet=TRUE), silent=TRUE)
        t1 <- proc.time()[3]
        RUNNINGTIMES[k, 2] <- t1 - t0
        if (class(res2) == "try-error") res2 <- matrix(Inf, 4, n - 1)
        t0 <- proc.time()[3]
        ## parametric boostrap
        res3 <- try(chronosCI(chr, op.out, 100, 3, cal, quiet=TRUE), silent=TRUE)
        t1 <- proc.time()[3]
        RUNNINGTIMES[k, 3] <- t1 - t0
        if (class(res3) == "try-error") res3 <- matrix(Inf, 4, n - 1)
        BT <- branching.times(drop.tip(tr, "outgroup"))
        OUT[[k]] <- cbind(BT, t(res1), t(res2), t(res3))
        ## save results at each loop
        save(OUT, RUNNINGTIMES, file = OUTFILE)
    }
}
