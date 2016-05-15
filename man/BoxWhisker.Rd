\name{boxWhisker}
\alias{boxWhisker}
\title{A function to create a box and whisker plot by group ID
}
\description{
A function to create a box and whisker plot by group ID.
}
\usage{
boxWhisker(data, groups, palette = "Paired")
}
\arguments{
  \item{data}{
a (non-empty) matrix of data values
}
  \item{groups}{
a (non-empty) vector of group IDs with length equal to the number of rows in data
}
  \item{palette}{
A color palette for plotting.  The default is 'Paired.'  See colorbrewer2.org for alternatives.
}
}
\value{
Returns a box-whisker plot of the data by group ID.
}
\examples{
data(Nuclei)
data(Groups)
boxWhisker(Nuclei, Groups)

#changing the color palette

boxWhisker(data = Nuclei, groups = Groups, palette = 'Set1')
}

