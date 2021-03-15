#########################
### Quantile Function ###
#########################

qage <- function(p=0.5, ages,  method="StraussSadler", k=min(length(ages),5)) {
 	
	if (length(ages) < 2) {
		stop("More than one fossil age is needed for inference");
	}
 	
 	ages <- sort(ages, decreasing=TRUE);
	
	if (method=="StraussSadler") {
		n <- length(ages);
		range <- max(ages)-min(ages);
		AGE <- range*(1-p)^-(1/n) + min(ages);
	} else if (method=="Beta") {
		n <- length(ages);
		range <- max(ages) - min(ages);
		X <- range/qbeta(p=p, shape1=n, shape2=1, lower.tail=FALSE);
		AGE <- X + min(ages);
	} else if (method=="Solow") {
		AGE <- ages[1] + (ages[1]-ages[2])*p/(1-p);
	} else if (method=="NorrisPenG") {
		PenG <- ages[1] - ages[2];
		UltG <- exp(qlogis(p=p, location=log(PenG)));
		AGE <- UltG + ages[1];
	} else if (method=="NorrisGLin") {
		GLin <- ages[1] - ages[2]
		UltG <- exp(qlogis(p=p, location=log(GLin/2)));
		AGE <- UltG + ages[1];
	} else if (method=="RobertsSolow") {
		# Select the k oldest observations (5 by default)
		ages <- ages[1:k];
		
		# Estimate the shape paremter of the joint Weibull distribution
		t <- numeric();
		
		for (i in 1:(k-2)) {
			t[i] <- (ages[1] - ages[k]) / (ages[1] - ages[i+1]);
		}
		
		v <- 1/(k-1) * sum(log(t));
		
		# Calculate Su (Robert & Solow 2003, = c(alpha) Solow 2005 eq. 18)
		Su <- (-log(1-p)/k)^-v;
		
		# But Su must be >= 1 (otherwise the estimate may be younger than older fossil) so:
		Su <- pmax(Su, 1); # use pmax instead of max to allow for vectors of p and Su
		
		# Calculate the quantile (Solow 2005 eq. 17) (Something seems wrong)	
		#AGE <- ( ages[1] - Su*ages[k]) / ( 1 - Su)
		
		# Calculate the quantile (Roberts & Solow 2003)
		AGE <- ages[1] + ( ages[1] - ages[k]) / ( Su - 1 );
	} else {
		stop("Method '", method, "' is not recognized.", call.=FALSE);
	}
	
	return(AGE);
}

