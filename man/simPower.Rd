\name{simPower}

\alias{simPower}

\title{
An internal function of binomPower, which actually calculates the p value
}

\description{
An internal function of binomPower, which actually calculates the p value.
}

\usage{
simPower(ndads,mm,vv,tau2,nperms,nbins)
}

\arguments{
  \item{ndads}{
a (non-empty) numeric value indicating the number of dads.
}
  \item{mm}{
a (non-empty) numeric value indicating the mean number of offspring per dad per bin (normal dist). mm must be less than vv. 
}
  \item{vv}{
a (non-empty) numeric value indicating the variance in offspring per dad per bin (normal dist). vv must be great than mm. 
}
  \item{tau2}{
a (non-empty) numeric value indicating the dad effect (narrow-sense heritability ~ tau2/(tau2+(pi/sqrt(3))^2)).
}
  \item{nperms}{
a (non-empty) numeric value indicating the number of bootstrap permutations to use for caluclating a p value.
}
  \item{nbins}{
a (non-empty) numeric value indicating the number of bins, data are pooled before analysis.
}
}


\value{
Returns a p value for a given set of conditions over a specificed number of bootstrap permutations.
}

\examples{
#not run
}
