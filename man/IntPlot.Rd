\name{IntPlot}
\alias{IntPlot}

\title{A function to visualize the results of a MANOVA
}
\description{
The function produces an interaction plot to demonstrate the results of a MANOVA using the function interaction.plot.  
}

\usage{
IntPlot(Scores, Cov.A, Cov.B, pvalues = rep(1, 8), int.pvalues = rep(1, 4))
}

\arguments{
  \item{Scores}{
A (non-empty) numeric matrix of principle component scores or raw data.
}
  \item{Cov.A}{
A (non-empty) bivariate factor vector indicating the factor for each row in Scores
}
  \item{Cov.B}{
A (non-empty) bivariate factor vector indicating the factor for each row in Scores
}
  \item{pvalues}{
An optional vector of p values for each covariate across Scores.  The length of pvalues must equal the number of columns in Scores times 2.
}
  \item{int.pvalues}{
An optional vector of p values for each interaction.  The length of int.pvalues must equal the number of columns in Scores.
}
}

\value{
Outputs a .pdf file of the interaction plot and key relating color and linetype to significance to the working directory.  
}

\seealso{
\code{\link{interaction.plot}}
}
\examples{
data(Scores)
data(CondA)
data(CondB)

pvals<-c(0.03,0.6,0.05,0.07,0.9,0.2,0.5,0.3)
int.pvals<-c(0.3,0.45,0.5,0.12)

IntPlot(Scores,CondA,CondB,pvalues=pvals, int.pvalues=int.pvals)
}

