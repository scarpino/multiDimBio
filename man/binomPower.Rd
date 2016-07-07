\name{binomPower}

\alias{binomPower}

\title{
Power analysis for estimating the heritability of a binomial trait
}

\description{
Performs a power analysis for estimating the heritability of a binomial trait.  This function can take a long time to run if either nsims or nperms is large.
}

\usage{
binomPower(ndads, mm, vv, tau2, nperms, nsims, nbins, alpha = 0.05, doPlot=FALSE)
}

\arguments{
  \item{ndads}{
a (non-empty) numeric value indicating the number of dads.
}
  \item{mm}{
a (non-empty) numeric value indicating the mean number of offspring per dad per bin (normal dist).  mm must be less than vv.
}
  \item{vv}{
a (non-empty) numeric value indicating the variance in offspring per dad per bin (normal dist). vv. must be greater than mm.
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
  \item{alpha}{
a (non-empty) numeric value indicating the cutoff for significant p values.
}
  \item{doPlot}{
a (non-empty) logical value indicating whether to plot the results of the power analysis.
}
}


\value{
Returns a list and an optional set of .pdfs (if doPlot==TRUE).  The list contains:
  \item{roc}{
a data.frame with the summarized results of the power analysis.
}

  \item{params}{
a numeric matrix with the paramater values.
}

  \item{results}{
a numeric matrix with the full results of the analysis.
}

}

\examples{
	ndads <- c(9,18)
	mm <- 4.629634
	vv <- 6.31339
	tau2 <- c(0,0.5)
	nperms <- 2
	nsims <- 2
	nbins <- 3
	doPlot <- TRUE
	binomPower(ndads,mm,vv,tau2,nperms,nsims,nbins,doPlot)
}
