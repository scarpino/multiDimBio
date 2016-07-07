\name{getP}

\alias{getP}

\title{
An internal function for getting empirical p values
}

\description{
Simulates p values.
}

\usage{
getP(ndads, mm, vv, tau2, nperms, nsims, nbins)
}

\arguments{
  \item{ndads}{
a (non-empty) numeric value indicating the number of dads.
}
  \item{mm}{
a (non-empty) numeric value indicating the mean number of offspring per dad per bin (normal dist). mm must be less than vv. 
}
  \item{vv}{
a (non-empty) numeric value indicating the variance in offspring per dad per bin (normal dist). vv must be greater than mm.
}
  \item{tau2}{
a (non-empty) numeric value indicating the dad effect (narrow-sense heritability ~ tau2/(tau2+(pi/sqrt(3))^2)).
}
  \item{nperms}{
a (non-empty) numeric value indicating the number of bootstrap permutations to use for caluclating a p value.
}
  \item{nsims}{
a (non-empty) numeric value indicating the number of simulations to run per parameter combination.
}
  \item{nbins}{
a (non-empty) numeric value indicating the number of bins, data are pooled before analysis.
}
}


\value{
Returns a vector of simulated p values.  The list contains:

}

\examples{
	ndads <- c(9,18)
	mm <- 4.629634
	vv <- 6.31339
	tau2 <- c(0,0.5)
	nperms <- 2
	nsims <- 2
	nbins <- 3
	getP(ndads, mm, vv, tau2, nperms, nsims, nbins)
}
