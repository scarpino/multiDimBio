\name{partialF}
\alias{partialF}

\title{
	A function to compute partial F statistics
}
\description{
	This is an internal function used in FSelect.  It can only be used for two groups.  The partial F statistic is the additional contribution to the model from adding one more trait.
}
\usage{
partialF(m.lda, GROUP, T_pm1)
}

\arguments{
  \item{m.lda}{
  	An object of class 'lda'
}
  \item{GROUP}{
  	A factor vector indicating group membership
}
  \item{T_pm1}{
  	The F statistic calculated for a discriminant analysis with only the most informative trait.
  	}
}

\value{
Returns a partial F statistic
}

\references{
	Habbema J, Hermans J (1977). Selection of Variables in Discriminant Analysis by F-Statistics and Error Rate. Technometrics, 19(4), 487 - 493.
}

\seealso{
\code{\link{FSelect}}
}
\examples{
#Internal function used in FSelect

data(Nuclei)
data(Groups)

NPC<-floor(ncol(Nuclei)/5)

DAT.comp<-completeData(Nuclei, n_pcs = NPC) 
Groups.use<-c(1,2) 
use.DAT<-which(Groups==Groups.use[1]|Groups==Groups.use[2])

DAT.use<-Nuclei[use.DAT,]
GR.use<-Groups[use.DAT]

traitA<-2

mlda<-MASS::lda(GR.use~DAT.use[,traitA])

F1<-partialF(mlda,GR.use,0)

traitB<-1

mlda2<-MASS::lda(GR.use~DAT.use[,c(traitA,traitB)])

partialF(mlda2,GR.use,F1)

}
