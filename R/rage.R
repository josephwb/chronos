##################################
### Random Clade Age Generator ###
##################################

# This function uses the quantile function above to generate random ages from the corresponding models.

rage <- function(ages, n, ...) {
	ages <- as.matrix(ages);

	if (ncol(ages) > 2) {
		stop("'ages' must be a vector or a two-column matrix");
	}

	# When fossil ages are known exactly
	if (ncol(ages) == 1) {
		rA <- qage(p=runif(n), ages=ages[,1], ...);
	} else if (ncol(ages) == 2) {
		# When fossil ages are known from time intervals (max-min)

		# Check that bounds are correct for all fossils
		if (any(apply(ages, 1, diff) < 0)) {
			stop("Detected at least one fossil with older bound < younger bound");
		}
	
		rA <- numeric(length=n);
		for(i in 1:n) {
			A <- runif(nrow(ages), min=ages[,1], max=ages[,2]);
			rA[i] <- qage(p=runif(1), ages=A, ...);
		}
	}
	return(rA);
}

