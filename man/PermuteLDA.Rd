\name{PermuteLDA}
\alias{PermuteLDA}

\title{
	A function to determine whether two groups are in statistically different locations in multivariate space See Collyer and Adams 2007
}

\description{
The function calculates the multivariate distance between two groups across all traits and determines whether they differ signifcantly using a Monte Carlo randomization test.  The Monte Carlo randomization creates a null distribution by randomizing the residual deviation from the group mean across all individuals.  This method controls for heteroscedasticity and was designed by Collyer and Adams (2007) for use in analyzing data sets that have sparse groups sizes relative to the number of traits.
}

\usage{
PermuteLDA(Data, Groups, NPerm, Missing.Data = "Complete")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
  	A (non-empty), numeric matrix of data values
}
  \item{Groups}{
	A (non-empty), vector indicating group membership.  
}
  \item{NPerm}{
  	The number of permutations used to generate the null distribution.  The default is 100.
}
   \item{Missing.Data}{
The method used to handle missing data.  The default, 'Complete' will use CompleteData to impute missing data, setting Missing.Data='Remove' will remove all individuals with missing data.  FSelect cannot handle missing data.
}
}

\details{
Determining the statistical significance of a discriminate function analysis along with performing that analysis on sparse data sets, e.g. many traits observed on comparatively few individuals, is a challenge.	Collyer and Adams (2007) developed a Monte Carlo based algorithm for addressing both of those issues. Briefly, the test uses the underlying Var/Cov structure of the data and randomizes the group membership to calculate a null distribution. This test simultaneously controls for heteroscedasticity, a common problem in sparse data sets and allows the approximation of a p-value for the test. For the original implementation and formulation of the method see Collyer and Adams (2007) or http://www.public.iastate. edu/~dcadams/software.html. Unlike the FSelect implementation, PermuteLDA will work properly with an arbitrary number of groups. The time required to run the algorithm is non-linear in the number of groups.
}

\value{
Returns a data frame with four columns and the number of groups choose 2 rows.  Each row is a pairwise comparison between groups.  The column 'Pr' is the p value to reject the null hypothesis of no difference (a value in 'Pr' < 0.05 would result in rejecting the hypothesis that the two groups are not different.  The column 'Distance' is the multivariate distance between the two groups.  
}
\references{
Collyer M, Adams D (2007). Analysis of Two - State Multivariate Phenotypic Change in Ecological Studies. Ecology, 88(3), 683 - 692.

For an implementation of the original method coded in R see http://www.public.iastate. edu/~dcadams/software.html.
}

\seealso{
\code{\link{PermuteLDA}}
}

\examples{
	data(Nuclei)
	data(Groups)
	PermuteLDA(Nuclei,Groups,50)
}
