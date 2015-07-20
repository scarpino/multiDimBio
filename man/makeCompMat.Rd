\name{makeCompMat}
\alias{makeCompMat}

\title{A function to create a pairwise comparison matrix
}
\description{
This function creates a pairwise comparison matrix for n groups.  All possible pairwise combinations are created, with rows in the matrix equal to the desired comparison.  
}

\usage{
makeCompMat(ng)
}

\arguments{
  \item{ng}{
A single number indicating the total number of unique groups
}
}

\value{
Returns a matrix with two columns and ng choose 2 rows.
}

\seealso{
\code{\link{PermuteLDA}}
}
\examples{
makeCompMat(3)

makeCompMat(4)

data(Groups)
NGroups<-length(unique(Groups))

makeCompMat(NGroups)
}

