\name{Power}
\alias{Power}

\title{
	A function to estimate the error rate for FSelect and PermuteLDA.
}
\description{
Methods are implemented to compute the statistical power, in terms of the type II error rate, based on anticipated sample and effect sizes for FSelect() and PermuteLDA(). By default the power of both tests are determined by iterating over a range of effect and sample sizes. The default settings were selected to be representative of many behavioral genetic studies; however, users can input alternative sample and effect sizes.   For high values of trials this function can be very slow.
}
\usage{
Power(func = "PermuteLDA", N = "DEFAULT.N", effect.size = "DEFAULT.e", trials = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{func}{
  	A character string indicating which function to compute the power for, can be either 'PermuteLDA' or 'FSelect'
}
  \item{N}{
  	A (non-empty) vector of group sizes.  The lenght of N must be greater than 1 and tha minimum group size for 'FSelect' can not be less than 6.  The size of each group is N/2.
}
  \item{effect.size}{
  	A (non-empty) vector or single value of effect sizes.
}
  \item{trials}{
  	A number indicating the number of trials for each combination of N and effect.size to calculate the power.
}
}
\details{
The algorithm for the power analysis proceeds as follows:
1. Input sample and effect sizes2. Set the number of significant effects, e to 0. Note - Total number of traits is fixed at 63. Draw random deviates for the given sample size for 6 traits. Note - All traits not significant under this iteration are drawn from a N(0,1) distribution.4. Perform either FSelect() or PermuteLDA() and record the results. 5. Return to step 3 N times, recording the results each time. Note - N is set using thetrials input6. If e<5 return to step 2 and set the number of significant effects to e+17. Proceed to the next combination of sample and effect size.8. Output the results for each combination of sample and effect size as a function of the number of significant traits.
}
\value{
Outputs a .pdf file for each effect size.
}
\seealso{
\code{\link{PermuteLDA}},\code{\link{FSelect}}
}
\examples{
#not run
#Power(func = 'FSelect', N=c(6,8), effect.size=0.5, trials = 2)
}
