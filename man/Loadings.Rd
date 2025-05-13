\name{Loadings}

\alias{Loadings}

\title{
	A function to visualize trait loadings onto discriminant function and principle component axes
}
\description{
This function produces barplots representative of the contribution of a particular trait or variable to either a discriminant function or principle component axis. 
}
\usage{
Loadings(DATA, GROUPS, method = c("PCA", "LDA"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DATA}{
  	A (non-empty) numeric matrix with trait values
}
  \item{GROUPS}{
  	A (non-empty)factor vector indicating the group membership of each row in DATA
}
  \item{method}{
  	An optional list indicating whether the results for a principle component analysis, 'PCA', or linear discriminant analysis, 'LDA' should be performed.
}
}

\value{
Outputs a list with values and plots for each test listed in method.
}

\seealso{
\code{\link[pcaMethods:pcaMethods]{pca}}, \code{\link[MASS]{lda}}
}
\examples{
data(Nuclei) 
data(Groups) 
Loadings(Nuclei, Groups, method=c("PCA", "LDA"))
}

