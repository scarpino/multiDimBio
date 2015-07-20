\name{h2Estimate}

\alias{h2Estimate}

\title{
Estimates the heritability of a binomial trait
}

\description{
Estimates the narrow-sense heritability of a binomial trait and calculates a p value by randomization.
}

\usage{
h2Estimate(data,nreps=1000)
}

\arguments{
  \item{data}{
a (non-empty) numeric matrix with three columns.  The first two should contain the trait data (number of occurances of each outcome type) and the third should contain the group ids.
}
  \item{nreps}{
a (non-empty) numeric value indicating the number of resamples to perform when calculating the emperical p value.
}
}

\details{
Estimates the narrow-sense heritability of a binomial trait.  This function works by fitting two models, one with and one without a random-effect of sire.  These models are compared by randomizing the sire ids nreps times and re-fitting the model.  For each of the nreps model pairs, a deviance is calculated and a "p value" estimated by comparing that distribution of deviance to the observed.  The heritability is approximatly tau2/(tau2+(pi/sqrt(3))^2), where tau2 is the random-effect variance due to sire.
}

\value{
Returns a list.  The list contains:
  \item{h2}{
The estimated narrow-sense heritability.  The narrow-sense heritability is approximatly tau2/(tau2+(pi/sqrt(3))^2), where tau2 is the random-effect variance due to sire.
}

  \item{pval}{
The probability that the best-fit model includes an extra variance term for sire (random effect of dad).  The value is calculated by comparing the deviances from nreps number of randomized model comparisions. 
}

  \item{deviance}{
The deviance between a null model without a random effect of dad and a model with.
}

  \item{sim}{
The simulated deviances used in calculating the p value in pval.
}

  \item{obsMod}{
The glmer model object resulting from the observed data.
}

}

\examples{
	ndads <- 18
	mm <- 4
	vv <-  6
	tau2 <- 1.5
	nbins <- 3
	
	mylogit <- function(x) log(x/{1-x})
	ilogit <- function(x) 1/{1+exp(-x)}
	
	swimprob <- ilogit(rnorm(ndads, 0, sqrt(tau2)))
	mytable <- NULL
	for(i in 1:ndads) {
		bincounts <- pmax(1,rnbinom(nbins, mu = mm, size = mm^2/{vv-mm}))
		swim <- rbinom(3, bincounts,swimprob[i])
		set <- bincounts - swim
		theserows <- data.frame(set=set,swim=swim, Dad = i, Bin = 1:nbins)
		mytable <- rbind(mytable, theserows)
	}
	
	est <- h2Estimate(mytable,nreps=10)
	
	print(est$h2)
}
