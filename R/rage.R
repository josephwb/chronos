##################################
### Random Clade Age Generator ###
##################################

# This function uses the quantile function above to generate random ages from the
# corresponding models.

##' @title Random clade age generator
##' @description This function uses the quantile function above to generate
##' random clade ages from the corresponding models.
##' @param n Number of samples to be drawn. Default n = 1000.
##' @param ages Either a vector of fossil ages or a matrix with two columns: the
##' first with the minimum age bounds (upper stratigraphic bounds) and the
##' second with the maximum age bounds (lower stratigraphic bounds) of each
##' fossil (in rows). A minimum of 2 are required.
##' @param ... Other options passed to \code{\link{qage}}.
##' @details If ages are known exactly, only min.ages is used. If some or all
##' ages have uncertainty, typically upper and lower bounds defined by
##' bracketing geological strata, age maximums and minimums are set for each
##' fossil.
##' @return A numeric vector of length \code{n} representing simulated clade
##' ages.
##' @importFrom stats runif
##' @examples
##' \dontrun{
##'   # Example using the default method of Strauss & Sadler (1989)
##'   hist(rage(n=10000, ages=c(50, 30, 25, 14, 3.5)), breaks=500, freq=FALSE,
##'     xlim=c(50,120), xlab="Clade Age",
##'     main="Strauss & Sadler Example\nn=1000, ages=c(50, 30, 25, 14, 3.5)");
##'   graphics::curve(dexp(x-50, rate=0.07), col='blue', lwd=2, add=TRUE);
##'   graphics::legend("right", legend="dexp", col="blue", lty=1, lwd=2, box.lty=0);
##'   
##'   # Example using the method of Solow (2003)
##'   hist(rage(n=10000, ages=c(50, 30, 25, 14, 3.5), method="Solow"),
##'     breaks=500, freq=FALSE, xlim=c(50,120), xlab="Clade Age",
##'     main="Solow Example\nn=1000, min.ages=c(50, 30, 25, 14, 3.5)");
##'   graphics::curve(dexp(x-50, rate=0.04), col='blue', lwd=2, add=TRUE); # heavy on intermediate values
##'   graphics::curve(dlnorm(x-50, meanlog=log(50-30), sdlog=pi/sqrt(3)), col='orange', lwd=2, add=TRUE);
##'   graphics::legend("right", legend=c("dexp", "dlnorm"), col=c("blue", "orange"), lty=1, lwd=2, box.lty=0);
##'   
##'   # Strauss & Sadler (1989) method incorporating fossil age uncertainty
##'   hist(rage(n=10000, ages=cbind(c(50, 30, 25, 14, 3.5), c(56, 35, 25, 14, 6))),
##'     breaks=500, freq=FALSE, xlim=c(50,120), xlab="Age",
##'     main="Solow Example (with fossil age uncertainty)n=1000,\nages=cbind(c(50, 30, 25, 14, 3.5), c(56, 35, 25, 14, 6))")
##'   graphics::curve(dlnorm(x-50, meanlog=2.6, sdlog=0.9), col='blue', lwd=2, add=TRUE);
##'   graphics::legend("right", legend="dlnorm", col="blue", lty=1, lwd=2, box.lty=0);
##'   }
##' @export
rage <- function(ages, n=1000, ...) {
	ages <- as.matrix(ages);

	if (ncol(ages) > 2) {
		stop("'ages' must be a vector or a two-column matrix.");
	}
	if (n < 0) {
	  stop("n must be a positive integer.")
	}

	# When fossil ages are known exactly
	if (ncol(ages) == 1) {
		rA <- qage(p=stats::runif(n), ages=ages[,1], ...);
	} else if (ncol(ages) == 2) {
		# When fossil ages are known from time intervals (max-min)

		# Check that bounds are correct for all fossils
		if (any(apply(ages, 1, diff) < 0)) {
			stop("Detected at least one fossil with older bound < younger bound.");
		}
	
		rA <- numeric(length=n);
		for(i in 1:n) {
			A <- stats::runif(nrow(ages), min=ages[,1], max=ages[,2]);
			rA[i] <- qage(p=stats::runif(1), ages=A, ...);
		}
	}
	return(rA);
}

