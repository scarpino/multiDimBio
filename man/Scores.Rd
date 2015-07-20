\name{Scores}
\alias{Scores}
\docType{data}
\title{Principle component scores based on the data in Nuclei
}
\description{
Principle component scores were computed using PPCA for the data set Nuclei.  
}
\usage{
data(Scores)
}

\format{
  A numeric matrix with 4 columns and the same number of rows as Nuclei.  There are no missing values.
}


\source{
The data are provided courtesy of David Crews at the University of Texas at Austin.
}
\references{
Crews, D, R Gillette, SV Scarpino, M Manikkam, MI Savenkova, MK Skinner. 2012.
Epigenetic Transgenerational Alterations to Stress Response in Brain Gene Networks and 
Behavior. Proc. Natl. Acad. Sci. USA. 109 (23). 9143 - 9148.
}

\examples{
data(Scores)

data(Nuclei)

SCORES<-PPCA(Nuclei)@scores
}
