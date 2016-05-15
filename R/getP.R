getP <- function(ndads, mm, vv, tau2, nperms, nsims, nbins) {
	ps <- c()
    for (i in 1:nsims) {
		p.i <- simPower(ndads, mm, vv, tau2, nperms, nbins)
		ps <- c(ps, p.i)
    }
	return(ps)
}  #end getP