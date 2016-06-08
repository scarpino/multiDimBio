#' An internal function for getting empirical p values
#' 
#' Simulates p values.
#' 
#' 
#' @param ndads a (non-empty) numeric value indicating the number of dads.
#' @param mm a (non-empty) numeric value indicating the mean number of
#' offspring per dad per bin (normal dist).
#' @param vv a (non-empty) numeric value indicating the variance in offspring
#' per dad per bin (normal dist).
#' @param tau2 a (non-empty) numeric value indicating the dad effect
#' (narrow-sense heritability ~ tau2/(tau2+(pi/sqrt(3))^2)).
#' @param nperms a (non-empty) numeric value indicating the number of bootstrap
#' permutations to use for caluclating a p value.
#' @param nsims a (non-empty) numeric value indicating the number of
#' simulations to run per parameter combination.
#' @param nbins a (non-empty) numeric value indicating the number of bins, data
#' are pooled before analysis.
#' @return Returns a vector of simulated p values.  The list contains:
#' @examples
#' 
#' 	ndads <- c(9,18)
#' 	mm <- 4.629634
#' 	vv <- 6.31339
#' 	tau2 <- c(0,0.5)
#' 	nperms <- 2
#' 	nsims <- 2
#' 	nbins <- 3
#' 	getP(ndads = ndads, mm = mm, vv = vv, tau2 = tau2, nperms = nperms, nsims = nsims, nbins = nbins)
#' 
#' @export getP
getP <- function(ndads, mm, vv, tau2, nperms, nsims, nbins) {
	ps <- c()
    for (i in 1:nsims) {
		p.i <- simPower(ndads, mm, vv, tau2, nperms, nbins)
		ps <- c(ps, p.i)
    }
	return(ps)
}  #end getP
