\name{ldaPlot}

\alias{ldaPlot}

\title{
	A function to visualize the results of a discriminant analysis
}
\description{
The function takes as input the traits and group IDs and will perform a discriminate function analysis and visualize the results. For the pair-wise comparison of groups we use density histograms with points along the x-axis denoting the actual data, Figure 3 For multi-group comparisons we plot a bivariate scatter for all pairwise combinations of discriminate axes. The color of plotting symbols can be altered using the palette argument and the axes comparisons (with max n = number of groups - 1).
}

\usage{
ldaPlot(Data, Groups, palette = "BrBG", axes = c(1, 2, 2, 3, 1, 3))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
  	A (non-empty), numeric matrix of data values
}
  \item{Groups}{
	A (non-empty), vector indicating group membership.  Length(unique(Group))==2
}
 	\item{palette}{
	A color palette for plotting.  The default is 'Paired.'  See colorbrewer2.org for alternatives.
}
  \item{axes}{
  	A numeric vector describing which axes to compare.  For example, axes=c(1,2) will on produce a single plot comparing the first and second axis.
  	}
}

\value{
Outputs a .pdf of the results
}


\seealso{
\code{\link{lda}}
}

\examples{
data(Nuclei)
data(Groups) 
ldaPlot(Nuclei, Groups, palette='BrBG', axes=c(1,2,2,3,1,3))
}
