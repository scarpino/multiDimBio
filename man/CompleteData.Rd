\name{completeData}

\alias{completeData}

\title{
Function to impute missing data.
}

\description{
This function imputes missing data using a probabilistic principle component analysis framework and is a wrapper around functions implemented in the pcaMethods package (Stacklies et al. 2007), was proposed by Troyanskaya et al 2001 and is based on methods developed in Roweis 1997.
}

\usage{
completeData(data, n_pcs, cut.trait = 0.5, cut.ind = 0.5, show.test = TRUE)
}

\arguments{
  \item{data}{
a (non-empty) numeric matrix of data values.
}
  \item{n_pcs}{
a (non-empty) numeric value indicating the desired number of principle component axes.
}
  \item{cut.trait}{
a number indicating the maximum proportion of missing traits before an individual is removed from data.  A value of 1 will not remove any individuals and 0 will remove them all.
}
  \item{cut.ind}{
a number indicating the maximum proportion of individuals missing a trait score before that trait is removed from data.  A value of 1 will not remove any traits and 0 will remove them all.
}
  \item{show.test}{
a logical statement indicating whether a diagnostic plot of the data imputation should be returned.
}
}


\value{
Returns a list with two entries. 
 \item{complete_dat}{
an object of class matrix with missing values imputed using a probabilistic principle component framework.
}
\item{plots}{
a list of plots stored as grid plots.
}

}

\references{
Roweis S (1997). EM algorithms for PCA and sensible PCA. Neural Inf. Proc. Syst., 10, 626 - 632.

Stacklies W, Redestig H, Scholz M, Walther D, Selbig J (2007). pcaMethods - a Bioconductor package providing PCA methods for incomplete data. Bioinformatics, 23, 1164 - 1167.

Troyanskaya O, Cantor M, Sherlock G, Brown P, Hastie T, Tibshirani R, Botstein D, Altman R (2001). Missing value estimation methods for DNA microarrays. Bioinformatics, 17(6), 520 - 5252.
}

\seealso{
\code{\link[pcaMethods:pcaMethods]{pcaMethods}}, \code{\link[pcaMethods:pcaMethods]{pca}}
}

\examples{
data(Nuclei)
npcs<-floor(ncol(Nuclei)/5)

length(which(is.na(Nuclei))==TRUE)

dat.comp<-completeData(data = Nuclei, n_pcs = npcs)

length(which(is.na(dat.comp))==TRUE)
}
