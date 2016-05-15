\name{FSelect}

\alias{FSelect}

\title{A Function to perform step-wise discriminant analysis using F statistics}

\description{
Select data using a F tests
}

\usage{
FSelect(Data, Group, target, p.adj.method = "holm", Missing.Data = "Complete")
}

\arguments{
  \item{Data}{
  	A (non-empty), numeric matrix of data values
}
  \item{Group}{
	A (non-empty), vector indicating group membership.  Length(unique(Group))==2
}
  \item{target}{
The number of desired traits. Target cannot be greater than the number of columns in Data
}
  \item{p.adj.method}{
The method used to control for false discovery.  The default setting is 'holm'
}
  \item{Missing.Data}{
The method used to handle missing data.  The default, 'Complete' will use completeData to impute missing data, setting Missing.Data='Remove' will remove all individuals with missing data.  FSelect cannot handle missing data.
}
}

\value{
FSelect saves the details of the final lda model as an object in the working directory.

FSelect returns list containing at least the following components:

\item{Selected}{
	An ordered list indicating which columns were selected.	
}
\item{F.Selected}{
	An ordered list containing the F statistics for each column indicated in Selected.  
}
\item{PrF}{
	An ordered list containing the p values for each column indicated in Selected.
}
\item{PrNotes}{
	A string indicating which method was used to control for multiple comparisons
}

}

\references{
Costanza M, Afifi A (1979). Comparison of Stopping Rules in Forward Stepwise Discriminant Analysis. Journal of the American Statistical Association, pp. 777 - 78

Habbema J, Hermans J (1977). Selection of Variables in Discriminant Analysis by F - Statistics and Error Rate. Technometrics, 19(4), 487 - 493.

Jennrich R (1977). Stepwise discriminant analysis, volume 3. New York Wiley Sons.
}

\seealso{
\code{\link{completeData}}
}

\examples{
data(Nuclei)
data(Groups)
npcs<-floor(ncol(Nuclei)/5)

dat.comp <- completeData(data = Nuclei, n_pcs = npcs) 
groups.use <- c(1,2) 
use.dat <- which(Groups==groups.use[1]|Groups==groups.use[2])

dat.use <- Nuclei[use.dat,]
GR.use <- Groups[use.dat]

#not run
#FSelect(DAT.use,GR.use,3)
}
