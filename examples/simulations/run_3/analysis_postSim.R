fls <- unzip("chronos_simulations.zip")

O <- NULL
RT <- NULL
for (i in 1:54) {
    fl <- grep(paste0("row-", i, "(\\.rda|_merged\\.rda)"), fls, value = TRUE)
    load(fl)
    if (length(OUT) < 100) OUT[(length(OUT) + 1):100] <- OUT[length(OUT)]
    RT <- rbind(RT, RUNNINGTIMES)
    O <- c(O, OUT)
}

## get PARA below
para <- PARA[rep(1:nrow(PARA), each = 100), ]

unlink(fls) # clean-up files

######################################################################

RUNNINGTIMES <- RT
OUT <- O

NTIPS <- c(20, 50, 100) # number of species (taxa, sequences)
NCALPT <- c(1, 5, 10) # number of calibration points
S <- c(1e3, 1e4) # number of nucleotides in the sequences
MODEL <- c(1, 2, 4) # the models
mod4chronos <- c("clock", "correlated", "relaxed", "relaxed")
## create a data frame with the parameter combinations:
PARA <- expand.grid(MODEL = MODEL, NCALPT = NCALPT, S = S, NTIPS = NTIPS)

X <- NULL
for (i in seq_along(OUT)) {
    M <- OUT[[i]]
    calnodes <- which(apply(M, 1, var) < 1e-8)
    ## drop the calibration nodes
    M <- M[-calnodes, ]
    BT <- M[, 1L] # true branching times
    for (type in 1:3) {
        start <- 2 + (type - 1) * 4
        n50 <- sum(BT >= M[, start + 1L] & BT <= M[, start + 2L]) # 50%
        n95 <- sum(BT >= M[, start] & BT <= M[, start + 3L]) # 95%
        ## E: errors (take the midpoints of the 50% intervals)
        E <- abs(BT - (M[, start + 1L] + M[, start + 2L])/2)
        me <- mean(E, na.rm = TRUE) # mean error
        sde <- sd(E, na.rm = TRUE) # SD error
        ## mean interval widths:
        iw50 <- mean(abs(M[, start + 1L] - M[, start + 2L]))
        iw95 <- mean(abs(M[, start] - M[, start + 3L]))
        if (!is.finite(me)) me <- NA
        X <- rbind(X, c(n50, n95, me, sde, iw50, iw95))
    }
}

colnames(X) <- c("n50", "n95", "me", "sde", "iw50", "iw95")
## find the parameters
para <- PARA[rep(1:nrow(PARA), 100L), ] # NOT each = ...!
para <- para[rep(1:nrow(para), each = 3), ]
para$Bootstrap <- factor(rep(1:3, length.out = nrow(para)), labels = c("Nonparametric", "Semiparametric", "Parametric"))
para$Nfreenodes <- para$NTIPS - para$NCALPT
for (j in 1:4) para[[j]] <- factor(para[[j]])
names(para)[1:4] <- c("Model", "Ncal", "s", "n")
X <- cbind(para, X)
levels(X$Model) <- c("Clock", "Correlated", "Uncorrelated")
X$RT <- as.vector(t(RUNNINGTIMES)) # add running times

## violin plots
library(vioplot)

BOTTOMSTRIP <- function(h = NULL, basecolour = "yellow") {
    if (is.null(h)) h <- yinch((par("mai")[1] - 0.1) / 3)
    par(xpd = TRUE)
    on.exit(par(xpd = FALSE))
    basecolour <- col2rgb(basecolour)
    fac <- c(255, 230, 211)
    alpha <- 0.5 * 255
    cols <- character(3)
    for (i in 1:3) {
        x <- basecolour * fac[i] / 255
        cols[i] <- rgb(x[1], x[2], x[3], alpha, maxColorValue = 255)
    }
    ## s
    xl <- 0.5; xr <- xl + 3; yt <- par("usr")[3]; yb <- yt - h
    labs <- rep(expression(s == 10^3, s == 10^4), 9)
    for (i in 1:18) {
        rect(xl, yb, xr, yt, col = cols[1])
        text((xl + xr)/2, (yb + yt)/2, labs[i])
        xl <- xr
        xr <- xl + 3
    }
    ## Ncal
    xl <- 0.5; xr <- xl + 6; yt <- yb; yb <- yt - h
    labs <- rep(paste0("Ncal = ", c("1", "5", "10")), 3)
    for (i in 1:9) {
        rect(xl, yb, xr, yt, col = cols[2])
        text((xl + xr)/2, (yb + yt)/2, labs[i])
        xl <- xr
        xr <- xl + 6
    }
    ## n
    xl <- 0.5; xr <- xl + 18; yt <- yb; yb <- yt - h
    labs <- paste0("n = ", c("20", "50", "100"))
    for (i in 1:3) {
        rect(xl, yb, xr, yt, col = cols[3])
        text((xl + xr)/2, (yb + yt)/2, labs[i], font = 2)#, col = "white")
        xl <- xr
        xr <- xl + 18
    }
}

TOPLEGEND <- function() {
    psr <- par("usr")
    xx <- psr[2]/2
    yy <- psr[4] * (0.5 + 0.5/par("plt")[4])
    legend(xx, yy, legend = paste0(c("Nonp", "Semip", "P"), "arametric"),
           pch = 22, pt.bg = cols, xjust = 0.5, yjust = 0.5,
           horiz = TRUE, bty = "n", pt.cex = 2, xpd = TRUE)
}

LINESEPARATORS <- function()
    abline(v = seq(3.5, by = 3, length.out = 17), lwd = 0.5, lty = 3)

cols <- c("indianred", "lightgreen", "skyblue2")

MYPLOT <- function(formula, model = "Clock", ylab = "95% CI width", ...)
    vioplot(formula, subset(X, Model == model), las = 2, col = cols,
            xaxt = "n", xlab = "", ylab = ylab, xaxs = "i", yaxs = "i", ...)

pdf("sim3_iw95.pdf", 10, 10)
layout(matrix(1:3, 3))
for (i in 1:3) {
    MYPLOT(iw95 ~ Bootstrap + s + Ncal + n, levels(X$Model)[i])
    BOTTOMSTRIP()
    TOPLEGEND()
    LINESEPARATORS()
    mtext(paste0("(", letters[i], ")"), line = 2, at = -0.5, font = 2)
}
dev.off()

pdf("sim3_n95.pdf", 10, 10)
layout(matrix(1:3, 3))
for (i in 1:3) {
    MYPLOT(n95/Nfreenodes ~ Bootstrap + s + Ncal + n, levels(X$Model)[i],
           "Proportion of dates in the 95% CIs")
    BOTTOMSTRIP()
    TOPLEGEND()
    LINESEPARATORS()
    mtext(paste0("(", letters[i], ")"), line = 2, at = -0.5, font = 2)
}
dev.off()

pdf("sim3_RT.pdf", 10, 10)
layout(matrix(1:3, 3))
for (i in 1:3) {
    MYPLOT(log10(RT) ~ Bootstrap + s + Ncal + n, levels(X$Model)[i],
           "Running time (sec)",
           ##expression(log[10]*"(running time in sec)"),
           yaxt = "n")
    axis(2, at = 0:4, labels = c(1, 10, 100, 1000, "10,000"), las = 1)
    BOTTOMSTRIP()
    TOPLEGEND()
    LINESEPARATORS()
    mtext(paste0("(", letters[i], ")"), line = 2, at = -1.5, font = 2)
}
dev.off()


## lattice plots
library(latticeExtra)
library(grid)

xl <- "Proportion of true dates within inferred 95% CI"
yl <- "Percentage"
p1 <- useOuterStrips(histogram(~ n95/Nfreenodes | Bootstrap + Model, X, subset = Ncal == 1, xlab = xl, ylab = yl))
p10 <- useOuterStrips(histogram(~ n95/Nfreenodes | Bootstrap + Model, X, subset = Ncal == 10, xlab = xl, ylab = yl))

xl <- "Width of 95% CI"
pb1 <- useOuterStrips(histogram(~ iw95 | Bootstrap + Model, X, subset = Ncal == 1, xlab = xl, ylab = yl))
pb10 <- useOuterStrips(histogram(~ iw95 | Bootstrap + Model, X, subset = Ncal == 10, xlab = xl, ylab = yl))

x <- 0.025

##pdf("CI_prop.pdf", 8, 10)
print(p1, split = c(1, 1, 1, 2), more = TRUE)
print(p10, split = c(1, 2, 1, 2))
grid.text(paste0("(", letters[1:2], ")"), x, c(1, 0.5) - x, gp = gpar(font = 2))
##dev.off()

##pdf("CI_iw.pdf", 8, 10)
print(pb1, split = c(1, 1, 1, 2), more = TRUE)
print(pb10, split = c(1, 2, 1, 2))
grid.text(paste0("(", letters[1:2], ")"), x, c(1, 0.5) - x, gp = gpar(font = 2))
##dev.off()



### analysis of running times ###

AIC(lm(RT ~ (n + Model + s + Ncal + Bootstrap)^2, X)) # <- the best so far
AIC(lm(RT ~ (n + Model + Bootstrap)^2, X))
AIC(lm(RT ~ (I((as.numeric(as.character(n)))) + Model + s + Ncal + Bootstrap)^2, X)) # not far


foo <- function(x)
    round(c(Mean = mean(x), SD = sd(x), Min = min(x), Max = max(x)), 1)

library(xtable)

Z <- rbind(aggregate(X$RT, by = list(X$n), FUN = foo),
           aggregate(X$RT, by = list(X$Model), FUN = foo),
           aggregate(X$RT, by = list(X$Bootstrap), FUN = foo),
           aggregate(X$RT, by = list(X$Ncal), FUN = foo),
           aggregate(X$RT, by = list(X$s), FUN = foo))

caption <- "Running times (in seconds). For each grouping, the times were averaged over the other categories."

print(xtable(cbind(Z$Group, as.data.frame(Z$x)), digits = 1, caption = caption), include.rownames = FALSE, booktabs = TRUE, caption.placement = "top")

########################################################################
