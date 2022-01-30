## Objectives: assess model selection with PHI-IC and corTest (Tao et al.)
## !!! expected to take ~70 h for 100 replications (~42 min per replication)

source("../COMMON_CODE.R")
getIC <- function(x) attr(x, "PHIIC")$PHIIC
## library(pegas) # to check unique sequences below
source("../rate.CorrTest.R")

nrep <- 100L # number of replications
## rate.sim <- 1 # rate used by phangorn::simSeq(); finally set by simData() below
chrctr <- chronos.control(dual.iter.max = 100)
pmlctr <- pml.control(trace = 0L)

## best to not erase outputs from the previous simulations:
OUTFILE <- getOutfile("sim_ModelSelection")
cat("##", date(), "\n", file = OUTFILE) # creates the file
cat("model", "n", "s", "Ncalpt", "IC.clk", "IC.cor", "IC.rel", "rho\n",
    file = OUTFILE, append = TRUE)


## the parameters to assess
NTIPS <- c(20, 50, 100) # number of species (taxa, sequences)
NCALPT <- c(1, 5, 10)
S <- c(1e3, 1e4) # number of nucleotides in the sequences
MODEL <- c(1:2, 4) # seq_along(FUNTREESIM) # the studied models

## create a data frame with the parameters:
PARA <- expand.grid(MODEL = MODEL, NCALPT = NCALPT, S = S, NTIPS = NTIPS)

set.seed(1) # seed used for the simulations reported in the paper
## modify or delete the previous line to have different results

for (i in 1:nrow(PARA)) {
    n <- PARA$NTIPS[i]
    Ncalpt <- PARA$NCALPT[i]
    s <- PARA$S[i]
    model <- PARA$MODEL[i]

    np <- n + Ntip(outtree)

    for (j in 1:nrep) {
        tr <- simTreeWithOutgroup(n)
        trsim <- rescaleEdgeLength(FUNTREESIM[[model]](tr))
## graphically compare both trees:
## layout(matrix(1:2, 2))
## plot(tr); axisPhylo()
## plot(trsim); axisPhylo()

        cal <- getCalibrationPoints(Ncalpt)
        X <- simData(trsim, s, 1, .5)
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

        chr.clk <- chronos(tr4chr, model = "clock", calibration = cal, quiet = TRUE, control = chrctr)
        chr.cor <- chronos(tr4chr, calibration = cal, quiet = TRUE, control = chrctr)
        chr.rel <- chronos(tr4chr, model = "relaxed", calibration = cal, quiet = TRUE, control = chrctr)
        rho <- rate.CorrTest(tr.ml, "outgroup", 0, "titi")

        cat(model, n, s, Ncalpt, getIC(chr.clk), getIC(chr.cor), getIC(chr.rel),
            rho, file = OUTFILE, append = TRUE)
        cat("\n", file = OUTFILE, append = TRUE)
    }
}

cat("##", date(), "\n", file = OUTFILE, append = TRUE)

