\name{LandscapePlot}

\alias{LandscapePlot}

\title{
	A function to visualize the Functional Landscape of measured traits
}
\description{
This function plots a three-dimensional landscape of measured traits.  The peak heights are relative with respect to the input data. The width of each peak is controlled by the argument sigma and has only an aesthetic purpose. The 3D image is generated using the \code{\link[misc3d:drawScene]{drawScene}} and \code{\link[misc3d:surfaceTriangles]{surfaceTriangles}} .   
}
\usage{
LandscapePlot(Data, Groups=NULL, PDF=FALSE,LocPlot=FALSE,control=c(75,1,30))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
  	A (non-empty) numeric matrix with trait values
}
  \item{Groups}{
  	A (non-empty)factor vector indicating the group membership of each row in Data. If there is only a single group present in Data then Groups=NULL (default).
}
  \item{PDF}{
  	Logical controlling whether to output the results as a .pdf or a .jpeg.  The default (PDF=FALSE) will produce a .jpeg.  The file size for .pdf output can be large.
}
  \item{LocPlot}{
  	Logical controlling whether to output a .pdf naming the peaks according to the columns they represent.  The defaul is FALSE.
}
	\item{control}{
	An optional numeric vector setting the control parameters for persp.  control[1] = theta, control[2] = r, control[3] = phi
}
}

\value{
a list of plots stored as grid plots (or.pdf if PDF=TRUE) file for each column in data.
}

\examples{
data(Nuclei) 
data(Groups)
 
#plotting the first six columns
#not run
#LandscapePlot(Nuclei[,1:6], Groups=Groups)
}
