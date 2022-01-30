##########################################################################
## This R script runs the analyses of the Felidae data reported in:     ##
## "Confidence intervals in molecular dating by maximum likelihood"     ##
## by Paradis et al. (in revision for Methods in Ecology and Evolution) ##
##########################################################################

##########################################################################
## Minimum requirements:
##   ape (5.6-1)
##   phangorn (2.8)
##   chronos (0.1-1)
## (Note: current versions of ape and phangorn on CRAN are enough;
##  chronos should be installed from GitHub)
##   Connection from Internet to run ape::read.GenBank.
##   MUSCLE3 should be installed and working when called from ape::muscle
##     (see ?muscle in R).
##########################################################################

## The file 'nuccore_result.txt' is the result from a search on NCBI
## with the keywords "Felidae Johnson 2006 Science". This returned 998
## accessions which were saved in a file (menu "Send to:" at the bottom).

## read the accesion numbers:
x <- scan("nuccore_result.txt", what = "", sep = "\n")
accnum <- gsub(" .*$", "", x[c(FALSE, FALSE, TRUE)])

## load ape and download the sequences from GenBank:
library(ape)
X <- read.GenBank(accnum, chunk.size = 200)
## save the downloaded sequences to work later on them:
saveRDS(X, "X.rds")

## tree downloaded from ToL (not used in the paper)
tr <- read.tree("subtree-ottol-563159-Felidae.tre") # 93 tips
a <- stripLabel(tr$tip.label)
b <- unique(a)
todrop <- NULL
for (tax in b) {
    i <- which(a == tax)
    if (length(i) == 1) next
    todrop <- c(todrop, i[-1])
}
trb <- drop.tip(tr, todrop)
plot(trb)

## restore the data from GenBank:
X <- readRDS("X.rds")

## extract the gene (locus) name from each sequence:
GENE <- gsub(" gene,.*$", "", attr(X, "description"))
GENE <- gsub(" \\(.*\\)$", "", GENE)
GENE <- gsub("^.* ", "", GENE)
length(table(GENE)) # how many sequences for each gene

## we drop the data for the locus ZFX
z <- which(GENE == "ZFX")
SP <- attr(X, "species")
SP[z]
X <- X[-z]
SP <- SP[-z]
GENE <- GENE[-z]

## SP is the vector of species associated to each sequence:
table(SP)

## Ngenes: number of genes
Ngenes <- length(table(GENE))
uGENE <- unique(GENE)

names(X) <- SP

## create a list for each gene before alignment
XX <- vector("list", Ngenes)

## align each gene separately with MUSCLE 3 and store the alignments in XX
for (i in 1:Ngenes) XX[[i]] <- muscle(X[which(GENE == uGENE[i])], quiet = FALSE)

## save the aligned sequences
saveRDS(XX, "Xali.rds")

XX <- readRDS("Xali.rds") # restore the aligned sequences

## create the "supermatrix":
Y <- XX[[1]]
for (i in 2:Ngenes) Y <- cbind(Y, XX[[i]], fill.with.gaps = TRUE)

names(XX) <- uGENE

## show the supermatrix (better to save as PNG)
##pdf("alig.pdf", 16, 9)
png("alig.png", 1600*.9, 900*.9)
mar <- par("mar")
mar[2] <- 6
par(mar = mar, xpd = TRUE)
image(Y)
NCOL <- sapply(XX, ncol)
l <- 1; r <- NCOL[1]; b <- -3.6; t <- -1.6
COL <- rep(c("red", "green", "lightblue"), length.out = Ngenes)
for (i in 1:Ngenes) {
    rect(l, b, r, t, col = COL[i])
    text((r + l)/2, (b + t)/2, uGENE[i])
    l <- r + 1
    r <- r + NCOL[i + 1]
}
dev.off()

## see the base frequencies for each gene
dim(Y)
base.freq(Y, 1, 1)
sapply(XX, base.freq, freq = TRUE, all = TRUE)

## load phangorn, do NJ tree, define a simple ML model, and fit with GTR+G+I
## with stochastic tree search
library(phangorn)
phydat <- phyDat(Y)
tr <- NJ(dist.ml(phydat))
m <- pml(tr, phydat, k = 4) # define the ML phylogenetic model
m1 <- optim.pml(m, optNni = TRUE, optBf = TRUE, optQ = TRUE, optInv = TRUE,
                optGamma = TRUE, rearrangement = "stochastic")

## save the output from phangorn with the ML tree results
saveRDS(m1, "m1.rds")
m1 <- readRDS("m1.rds") # restore the output from phangorn

## get the data; same as: Y <- phyDat(Y)
Y <- m1$data

## get the outgroup:
outgroup <- labels(Y)[38:44]

## root the tree and drop the outgroup:
phy <- drop.tip(root(m1$tree, outgroup), outgroup)

## define the calibration point and the control paramters:
cal <- makeChronosCalib(phy, age.min = 16, age.max = 20)
ctr.chr <- chronos.control(dual.iter.max = 200)

## load the file from Tao et al's RateCorrTest
source("../sim/rate.CorrTest.R")
rate.CorrTest(root(m1$tree, outgroup), outgroup)

## fit the 3 main models in chronos():
chr <- chronos(phy, cal = cal, control = ctr.chr)
chr.rel <- chronos(phy, model = "relaxed", cal = cal, control = ctr.chr)
chr.clk <- chronos(phy, model = "clock", cal = cal, control = ctr.chr)

## load chronos
library(chronos)

## compute the CIs with the three bootstrap methods
ci.clk <- chronosCI(chr.clk, m1, cal = cal)
cib.clk <- chronosCI(chr.clk, m1, cal = cal, type = "nonparametric")
cic.clk <- chronosCI(chr.clk, m1, cal = cal, type = "parametric")

## combine them
CI <- list(ci.clk, cib.clk, cic.clk)

saveRDS(CI, "CI.rds") # save for the plot

CI <- readRDS("CI.rds") # restore the CIs

apply(sapply(CI, function(x) x[4, ] - x[1, ]), 2, summary) # 95% CIs
apply(sapply(CI, function(x) x[3, ] - x[2, ]), 2, summary) # 50% CIs

rates <- attr(chr, "rates")

## final chronogram with its CIs:
pdf("Felidae_3.pdf", 10, 20)
layout(matrix(1:3, 3))
par(mar = c(5, 4, 3, 0), xpd = TRUE, cex = .95)
for (i in 1:3) {
    o <- plot(chr.clk, edge.w = 1 * rates / max(rates), label.offset = .1, x.lim = 25)
    axisPhylo()
    text(o$x.lim[1], o$y.lim[2], paste0("(", letters[i], ")"), font = 2, cex = 1.2, adj = 2)
    mtext("Time (Ma)", 1, 3, at = 10)
    drawChronosCI(CI[[i]], legend = "bottomleft", h = .6)
}
dev.off()

######################################################################
## Bonus: rooted ML with phangorn

tr_rooted <- wpgma(dist.ml(phydat))

m_rooted <- pml(tr_rooted, phydat, k = 4) # define the ML phylogenetic model
m_rooted <- optim.pml(m_rooted, optNni = TRUE, optBf = TRUE, optQ = TRUE,
                optRooted=TRUE, optGamma = TRUE, optInv=TRUE,
                rearrangement = "stochastic")
plot(m_rooted)
axisPhylo()

saveRDS(m_rooted, "m_rooted.rds")

bs_rooted <- bootstrap.pml(m_rooted, multicore = TRUE, mc.cores = 4L, optRooted = TRUE)
saveRDS(bs_rooted, "bs_rooted.rds")

outgroup <- labels(Y)[38:44]
tree <- m_rooted$tree
tree <- drop.tip(tree, outgroup)
bs_trees <- lapply(bs_rooted, drop.tip, outgroup)
class(bs_trees) <- "multiPhylo"
tmp <- sapply(bs_trees, branching.times)
CI <- apply(tmp, 1, quantile, probs = c(0.025, 0.25, 0.75, 0.975))

plot(tree)
axisPhylo()
drawChronosCI(CI, h = .6)

## scaling with root age = 18 Ma
fac <- 18 / branching.times(tree)[1]
trb <- tree
trb$edge.length <- fac * tree$edge.length
CIb <- fac * CI

pdf("Felidae_rooted_ML.pdf", 10, 8)
par(xpd = TRUE)
plot(trb)
axisPhylo()
drawChronosCI(CIb, h = .6)
dev.off()
