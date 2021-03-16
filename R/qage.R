#########################
### Quantile Function ###
#########################

##' @title Quantile-based age estimation
##' @description Quantile function that computes the age corresponding to a
##' particular probability for the upper bound of a distribution of ages
##' (Strauss & Sadler 1989, Gingerich & Uhen 1998, Solow 2003). The method of
##' Strauss & Sadler (1989) assumes that the distribution of fossil ages is
##' uniform and their formula depends on the fossil ages range and the number of
##' fossil ages. The method of Solow (2003) is a general method for non-uniform
##' distributions and depends on the temporal gap between the oldest and the
##' second oldest fossil ages. Both methods assume that fossil ages are
##' independent samples from the same distribution (only relevant for the two
##' oldest ages for Solow's method), therefore, fossils should be as independent
##' as possible (ideally from different geological formations and different
##' regions).
##' 
##' In the particular case where there are only two fossil ages, Strauss &
##' Sadler's and Solow's methods converge to the same result; the quantile
##' functions are simply Xn/(1-P), and the likelihood function is 1/X.
##' @param p The desired probability level (0 < p < 1). A vector of
##' probabilities can be provided. Default p = 0.5
##' @param ages Either a vector of fossil ages or a matrix with two columns: the
##' first with the minimum age bounds (upper stratigraphic bounds) and the
##' second with the maximum age bounds (lower stratigraphic bounds) of each
##' fossil (in rows). A minimum of 2 are required
##' @param method The method for modelling age uncertainty. A number of options
##' is available:
##' \itemize{
##' \item{\code{"StraussSadler"}} {The (default) method of Strauss & Sadler 
##' (1989) assumes that the sampled fossil ages are uniformly distributed in 
##' time, and a warning is returned if a Kolmogorov-Smirnov test rejects the 
##' uniformity hypothesis.}
##' \item{\code{“Beta”}} {The method is a different parameterization of the
##' Strauss & Sadler method (Wang et al., 2009) that uses the qbeta function,
##' since the ratio between the observed maximum fossil age and the clade age
##' for the Strauss & Sadler model is distributed according to a Beta
##' distribution with parameters N and 1 (Wang et al. 2009); this should give
##' the same result as "StraussSadler".}
##' \item{\code{Solow}} {The method of Solow (2003) does not assume uniformity 
##' in sampled fossil ages, but is based on the two oldest ages only.}
##' \item{\code{NorrisPenG} or \code{NorrisGLin}} {The method of Norris et al.
##' (2015) based on the two oldest fossils only and the log-logistic
##' distribution. "NorrisPenG" is used when the precise phylogenetic placement
##' of fossils is not known, whereas "NorrisGLin" is used when one fossil from
##' each daughter lineage is used.}
##' \item{\code{"RobertsSolow"}} {The method of Roberts & Solow (2003) does not
##' assume uniformity in sampling and is based on the fact that the joint
##' distribution of the k oldest fossil ages can be modeled with the same
##' Weibull distribution.}
##' }
##' @param k The number of fossil ages to use in when
##' \code{method="RobertsSolow"}, in which only the k oldest fossils are used.
##' Default k = 5.
##' @return A numeric value (or vector of numeric values, if multiple p values
##' are provided) representing the age estimate of the clade origin given the
##' method a p value provided
##' @importFrom stats qbeta qlogis
##' @examples
##' \dontrun{
##'   # The following demonstrates how inferences depend on p and method
##'   qage(p=c(0.1, 0.5, 0.9), ages=c(54, 30, 25, 14, 5))
##'   qage(p=c(0.1, 0.5, 0.9), ages=c(54, 30, 25, 14, 5), method="Beta")
##'   qage(p=c(0.1, 0.5, 0.9), ages=c(54, 30, 25, 14, 5), method="Solow")
##'   }
##' @importFrom Rdpack reprompt
##' @references
##' \insertRef{Gingerich1998}{chronos}
##' 
##' \insertRef{Norris2015}{chronos}
##' 
##' \insertRef{Solow2003}{chronos}
##' 
##' \insertRef{Strauss1989}{chronos}
##' 
##' \insertRef{Wang2007}{chronos}
##' 
##' \insertRef{Wang2009}{chronos}
##' 
##' \insertRef{Wang2010}{chronos}
##' @export
qage <- function(p=0.5, ages,  method="StraussSadler", k=5) {
 	
	if (length(ages) < 2) {
		stop("More than one fossil age is needed for inference.", call.=FALSE)
	}
  if (p < 0 || p > 1.0) {
    stop("Valid values for p are: 0 < p < 1.")
  }
 	
 	ages <- sort(ages, decreasing=TRUE)
	
	if (method=="StraussSadler") {
		n <- length(ages)
		range <- max(ages) - min(ages)
		AGE <- range*(1-p)^-(1/n) + min(ages)
	} else if (method=="Beta") {
		n <- length(ages)
		range <- max(ages) - min(ages)
		X <- range/stats::qbeta(p=p, shape1=n, shape2=1, lower.tail=FALSE)
		AGE <- X + min(ages)
	} else if (method=="Solow") {
		AGE <- ages[1] + (ages[1]-ages[2])*p/(1-p)
	} else if (method=="NorrisPenG") {
		PenG <- ages[1] - ages[2]
		UltG <- exp(stats::qlogis(p=p, location=log(PenG)))
		AGE <- UltG + ages[1]
	} else if (method=="NorrisGLin") {
		GLin <- ages[1] - ages[2]
		UltG <- exp(stats::qlogis(p=p, location=log(GLin/2)))
		AGE <- UltG + ages[1]
	} else if (method=="RobertsSolow") {
		# Select the k oldest observations
	  if (k > length(ages)) {
	    print(paste0("NOTE: ", k, " oldest fossil specified, but only ",
        length(ages), " fossils available. Using all fossils in calculations."))
	    k <- length(ages)
	  }
	  ages <- ages[1:k]
		
		# Estimate the shape paremter of the joint Weibull distribution
		t <- numeric()
		
		for (i in 1:(k-2)) {
			t[i] <- (ages[1] - ages[k]) / (ages[1] - ages[i+1])
		}
		
		v <- 1/(k-1) * sum(log(t))
		
		# Calculate Su (Robert & Solow 2003, = c(alpha) Solow 2005 eq. 18)
		Su <- (-log(1-p)/k)^-v
		
		# But Su must be >= 1 (otherwise the estimate may be younger than older fossil) so:
		Su <- pmax(Su, 1); # use pmax instead of max to allow for vectors of p and Su
		
		# Calculate the quantile (Solow 2005 eq. 17) (Something seems wrong)	
		#AGE <- ( ages[1] - Su*ages[k]) / ( 1 - Su)
		
		# Calculate the quantile (Roberts & Solow 2003)
		AGE <- ages[1] + ( ages[1] - ages[k]) / ( Su - 1 )
	} else {
		stop("Method '", method, "' is not recognized.", call.=FALSE)
	}
	
	return(AGE)
}

