\name{ZTrans}
\alias{ZTrans}

\title{
	A function to convert data into a z-score
	}
\description{
This function converts the columns in a data matrix into z-scores.  The score is computed by subracting each observation in a column from the column mean and divding by the column standard deviation.  Each column is converted independently of the others missing values are ignored in the calculation.
}
\usage{

ZTrans(DATA)
}

\arguments{
  \item{DATA}{
  	A (non-empty) matrix with data values.  Columns should be different traits and rows unique observations of those traits
}
}

\value{
Returns a matrix with the same dimensions as DATA.
}

\seealso{
\code{\link{PercentMax}}, \code{\link{MeanCent}}
}
\examples{
data(Nuclei)

colMeans(Nuclei, na.rm=TRUE)

Nuclei.ZT<-ZTrans(Nuclei)

colMeans(Nuclei.ZT, na.rm=TRUE)
}

