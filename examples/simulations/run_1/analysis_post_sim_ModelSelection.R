## code to analyse the output, with figure production

OUTFILE <- "sim_ModelSelection010.out"
RES <- read.table(OUTFILE, TRUE)

## --> drop the results with lognormal uncorrelated model <-- ##
## (not reported in the paper)
RES <- RES[RES$model != 3, ]

selmod <- apply(RES[, c("IC.clk", "IC.cor", "IC.rel")], 1, which.min)
test <- ifelse(RES$rho > 0.5, "s", "ns")
finalsel <- ifelse(selmod == 1 & test == "s", 2, selmod)
table(RES$model, finalsel)
table(RES$model, test)
table(selmod, test)

for (n in NTIPS) {
    s <- RES$n == n
    print(table(RES$model[s], finalsel[s]))
}
for (nc in NCALPT) {
    s <- RES$Ncalpt == nc
    print(table(RES$model[s], finalsel[s]))
}
for (ss in S) {
    s <- RES$s == ss
    print(table(RES$model[s], finalsel[s]))
}

## x,y: coordinates of the topleft corner
RECT <- function(TAB, x, y, cex = 1, percentage = FALSE) {
    NR <- nrow(TAB); NC <- ncol(TAB)
    N <- rowSums(TAB)[1]
    Z <- sweep(TAB, 1, N, "/")
    col <- rgb(1 - Z, 1 - Z, 1)
    i <- 0L
    for (xx in x:(x + NC - 1L)) {
        for (yy in y:(y - NR + 1L)) {
            i <- i + 1L
            rect(xx, yy - 1, xx + 1, yy, col = col[i], border = NA)
            lab <- TAB[i]
            if (percentage) lab <- paste0(lab, "%")
            text(xx + 1/2, yy - 1/2, lab, col = ifelse(Z[i] <= .4, "black", "white"), cex = cex)
        }
    }
}

UPSTRIP <- function(x, y, label) {
    h <- 0.4
    rect(x, y, x + 3, y + h, col = "slategrey")
    text(x + 1.5, y + h/2, label, col = "white")
}

LOWSTRIP <- function(x, y, label) {
    h <- 0.4
    rect(x, y, x + 9, y - h, col = "slategrey")
    text(x + 4.5, y - h/2, label, col = "white")
}

LEFTSTRIP <- function(x, y, label) {
    w <- 0.4
    rect(x - w, y, x, y - 3, col = "slategrey")
    text(x - w/2, y  - 1.5, label, col = "white", srt = 90)
}


RES <- lapply(RES, function(x) if (is.integer(x)) factor(x))
finalsel <- factor(finalsel)

TAB <- table(RES$model, finalsel)
TAB <- round(100 * sweep(TAB, 1, rowSums(TAB), "/"), 1)

all(rowSums(TAB) == 100)
TAB[length(TAB)] <- 98.5
all(rowSums(TAB) == 100)

pdf("res_global_sim_modselec.pdf", 7.3, 5.7)
par(xpd = TRUE, mar = c(0, 4, 4, 0))
plot(NA, type = "n", asp = 1, xlim = c(0, 3), ylim = c(1, 4), axes = 0, ann = 0)
RECT(TAB, 0, 4, 1.5, TRUE)
rect(0, 1, 3, 4)
text(-0.5, 4:2 - 0.5, c("Strict clock", "Correlated\nlognormal", "Uncorrelated\ngamma"))
text(-.7, 4, "Simulated model", font = 2, adj = c(0.5, 1))
text(1:3 - 0.5, 4.2, c("Strict clock", "Correlated", "Uncorr. gamma"), xpd = TRUE)
text(1.5, 4.5, "Selected model", font = 2)
dev.off()

pdf("res_detailed_sim_modselec.pdf", 9, 5)
par(mar = rep(0.5, 4))
plot(NA, type = "n", asp = 1, xlim = c(0, 18), ylim = c(3, 12), axes = 0, ann = 0)
x <- 0; y <- 12
for (n in NTIPS) {
    for (Ncal in NCALPT) {
        sel <- RES$n == n & RES$Ncalpt == Ncal & RES$s == 1000
        TAB <- table(RES$model[sel], finalsel[sel])
        RECT(TAB, x, y)
        if (n == 20) LEFTSTRIP(x, y, paste("Ncal =", Ncal))
        if (Ncal == 1) UPSTRIP(x, y, paste("n =", n))
        y <- y - 3
    }
    y <- 12
    x <- x + 3
}
LOWSTRIP(0, 3, "s = 1000")
for (n in NTIPS) {
    for (Ncal in NCALPT) {
        sel <- RES$n == n & RES$Ncalpt == Ncal & RES$s == 1e4
        TAB <- table(RES$model[sel], finalsel[sel])
        RECT(TAB, x, y)
        if (Ncal == 1) UPSTRIP(x, y, paste("n =", n))
        y <- y - 3
    }
    y <- 12
    x <- x + 3
}
for (x in seq(3, 18, 3)) segments(x, 3, x, 12, lwd = 2)
for (y in c(6, 9)) segments(0, y, 18, y, lwd = 2)
LOWSTRIP(9, 3, "s = 10,000")
dev.off()


### same simulations but with milder parameters (and no lognormal) ###
OUTFILE <- "sim_ModelSelection011.out"
RES <- read.table(OUTFILE, TRUE)
## create TAB as above

pdf("res_global_sim_modselec_mildParameters.pdf", 5.7, 5)
par(xpd = TRUE, mar = c(0, 7.5, 4, 0))
plot(NA, type = "n", asp = 1, xlim = c(0, 3), ylim = c(1, 4), axes = 0, ann = 0)
RECT(TAB, 0, 4, 1.5, TRUE)
rect(0, 1, 3, 4)
text(-0.5, 4:2 - 0.5, c("Strict clock", "Correlated\nlognormal", "Uncorrelated\ngamma"))
text(-.7, 4, "Simulated model", font = 2, adj = c(0.5, 1))
text(1:3 - 0.5, 4.2, c("Strict clock", "Correlated", "Uncorr. gamma"), xpd = TRUE)
text(1.5, 4.5, "Selected model", font = 2)
dev.off()
