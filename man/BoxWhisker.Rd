\name{BoxWhisker}
\alias{BoxWhisker}
\title{A function to create a box and whisker plot by group ID
}
\description{
A function to create a box and whisker plot by group ID.
}
\usage{
BoxWhisker(DATA, GROUPS, palette = "Paired")
}
\arguments{
  \item{DATA}{
a (non-empty) matrix of data values
}
  \item{GROUPS}{
a (non-empty) vector of Group IDs with length equal to the number of rows in DATA
}
  \item{palette}{
A color palette for plotting.  The default is 'Paired.'  See colorbrewer2.org for alternatives.
}
}
\value{
Saves a box-whisker plot of the data by group ID as a .pdf in the working directory.
}
\examples{
data(Nuclei)
data(Groups)
BoxWhisker(Nuclei, Groups)

#changing the color palette

BoxWhisker(Nuclei, Groups, palette='Set1')
}

