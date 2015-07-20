\name{plotBinomPower}
\alias{plotBinomPower}
\title{A function to plot the results of a binomPower run
}
\description{
A function to plot the results of a binomPower run.
}
\usage{
plotBinomPower(datPlotBig,params,timeStamp=Sys.time())
}
\arguments{
  \item{datPlotBig}{
a (non-empty) matrix of data values, with columns trueTau, ndads, trueTau2
}
  \item{params}{
a (non-empty) matrix of parameter values, with columns mm and vv. 
}
  \item{timeStamp}{
A timeStamp.
}
}
\value{
Saves two plots of the binomPower analysis results a .pdf in the working directory.
}
\examples{
#not run
}

