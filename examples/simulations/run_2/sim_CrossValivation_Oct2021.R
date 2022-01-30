## Objectives: assess smoothing of the correlated model

source("../COMMON_CODE.R")

n <- 20
Ncalpt <- 10
s <- 1e4 # number of nucleotides in the sequences
rate.sim <- 1 # rate used by phangorn::simSeq()
chrctr <- chronos.control(dual.iter.max = 200)
pmlctr <- pml.control(trace = 0L)

set.seed(123)
tr <- simTreeWithOutgroup(n)
trsim <- rescaleEdgeLength(FUNTREESIM[[2]](tr))
cal <- getCalibrationPoints(10)
X <- simData(trsim, s, 1, 0.5)
utr <- unroot(trsim)
op.out <- optim.pml(pml(utr, X, rate = rate.sim), control = pmlctr)
tr.ml <- root(op.out$tree, "outgroup")
tr4chr <- drop.tip(tr.ml, "outgroup")

if (Ncalpt > 1) {
    for (i in 2:10) {
        nodetr <- cal$node[i]
        NODE <- getMRCA(tr4chr, PP[[nodetr - n]])
        if (NODE != nodetr) cal$node[i] <- NODE
    }
}

library(chronos)
cv.10 <- CV(tr4chr, 10^(-3:4), calibration = cal)
cal1 <- cal[1, ]
cv.1 <- CV(tr4chr, 10^(-3:4), calibration = cal1)

phiic <- numeric()
for (l in (lambda <- 10^(-3:4)))
    phiic <- c(phiic, attr(chronos(tr4chr, l, calibration = cal1, quiet = TRUE), "PHIIC")$PHIIC)
## => PHI-IC does not change a lot with lambda

xlab <- expression(lambda)
ylim <- range(c(cv.1[,2], cv.10[,2]))
xaxlabs <- expression(10^-3, 10^-2, 10^-1, 1, 10, 10^2, 10^3, 10^4)

pdf("CV_lambda.pdf", 8, 10)
layout(matrix(1:2, 2))
m <- par("mar")
m[3:4] <- 1
par(mar = m, las = 1)
plot(cv.1, NULL, "o", log = "xy", xlab = xlab, ylim = ylim, xaxt = "n")
points(cv.10, NULL, "o", col = "blue")
legend("topright", NULL, c(" 1 calibration point", "10 calibration points"),
       pch = 1, lty = 1, col = c("black", "blue"))
mtext("(a)", font = 2, at = 0.0001, cex = 1.2)
axis(1, at = lambda, labels = xaxlabs)
plot(lambda, phiic, "o", log = "x", xlab = xlab, ylab = expression(Phi*"IC"), xaxt = "n")
mtext("(b)", font = 2, at = 0.0001, cex = 1.2)
axis(1, at = lambda, labels = xaxlabs)
dev.off()

truetree <- drop.tip(tr, "outgroup")
chr1 <- chronos(tr4chr, 100, calibration = cal1, control = chrctr)
chr10 <- chronos(tr4chr, 10, calibration = cal, control = chrctr)
TR <- c(truetree, chr10, chr1)
names(TR) <- c("True tree", "10 cal., lambda = 10", "1 cal., lambda = 100")

chra <- chronos(tr4chr, calibration = cal1, control = chrctr)
chrb <- chronos(tr4chr, calibration = cal, control = chrctr)
TS <- c(truetree, chrb, chra)
names(TS) <- c("True tree", "10 cal., lambda = 1", "1 cal., lambda = 1")

yy <- seq(19, 16, length.out = 3)

pdf("LTT.pdf", 8, 10)
layout(matrix(1:2, 2))
mltt.plot(TR, dcol = FALSE, dlty = TRUE, legend = FALSE)
for (i in seq_along(yy)) {
    segments(-50, yy[i], -47.5, yy[i], lty = i)
    text(-47, yy[i], switch(i,"True tree",
                            expression("10 cal., "*lambda == 10),
                            expression("1 cal., "*lambda == 100)),
         adj = 0)
}
mtext("(a)", line = 1.5, font = 2, at = -57, cex = 1.2)
mltt.plot(TS, dcol = FALSE, dlty = TRUE, legend = FALSE)
for (i in seq_along(yy)) {
    segments(-50, yy[i], -47.5, yy[i], lty = i)
    text(-47, yy[i], switch(i,"True tree",
                            expression("10 cal., "*lambda == 1),
                            expression("1 cal., "*lambda == 1)),
         adj = 0)
}
mtext("(b)", line = 1.5, font = 2, at = -57, cex = 1.2)
dev.off()
## => more calibration points compensate for wrong smoothing parameter!
