load("sim_wrongmodel.rda")

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
para <- PARA[rep(1:nrow(PARA), each = nrep), ]
para <- para[rep(1:nrow(para), each = 3), ]
para$Bootstrap <- factor(rep(1:3, length.out = nrow(para)), labels = c("Nonparametric", "Semiparametric", "Parametric"))
para$Nfreenodes <- n - para$NCALPT

for (j in 1:3) para[[j]] <- factor(para[[j]])
names(para)[1:2] <- c("Model", "Ncal")
X <- cbind(para, X)
levels(X$Ncal) <- paste("Ncal =", levels(X$Ncal))
    ## levels(X$Model) <- c("Clock", "Correlated")

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
    ## Ncal
    xl <- 0.5; xr <- xl + 3; yt <- par("usr")[3]; yb <- yt - h
#    xl <- 0.5; xr <- xl + 6; yt <- yb; yb <- yt - h
    labs <- rep(paste0("Ncal = ", c("1", "5", "10")), 3)
    for (i in 1:3) {
        rect(xl, yb, xr, yt, col = cols[2])
        text((xl + xr)/2, (yb + yt)/2, labs[i])
        xl <- xr
        xr <- xl + 3
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
    abline(v = seq(3.5, by = 3, length.out = 2), lwd = 0.5, lty = 3)

MYPLOT <- function(formula, model = "Clock", ylab = "95% CI width", ...)
    vioplot(formula, subset(X, Model == model), las = 2, col = cols,
            xaxt = "n", xlab = "", ylab = ylab, xaxs = "i", yaxs = "i", ...)

cols <- c("indianred", "lightgreen", "skyblue2")

pdf("sim4.pdf", 12, 10)
layout(matrix(1:4, 2, 2, TRUE))
for (i in 1:2) {
    MYPLOT(iw95 ~ Bootstrap + Ncal, levels(X$Model)[i])
    BOTTOMSTRIP()
    TOPLEGEND()
    LINESEPARATORS()
    mtext(paste0("(", letters[i], ")"), line = 2, at = -0.1, font = 2)
}
for (i in 1:2) {
    MYPLOT(n95/Nfreenodes ~ Bootstrap + Ncal, levels(X$Model)[i],
           "Proportion of dates in the 95% CIs")
    BOTTOMSTRIP()
    TOPLEGEND()
    LINESEPARATORS()
    mtext(paste0("(", letters[i+2], ")"), line = 2, at = -0.1, font = 2)
}
dev.off()
