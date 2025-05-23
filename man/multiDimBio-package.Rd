\name{multiDimBio-package}
\alias{multiDimBio-package}
\alias{multiDimBio}
\docType{package}

\title{Multivariate Analysis and Visualization for Biological Data
}

\description{
Code to support a systems biology research program from inception through publication.  The methods focus on dimension reduction approaches to detect patterns in complex, multivariate experimental data and places an emphasis on informative visualizations.  The goal for this project is to create a package that will evolve over time, thereby remaining relevant and reflective of current methods and techniques.  As a result, we encourage suggested additions to the package, both methodological and graphical.
}

\details{
\tabular{ll}{
Package: \tab multiDimBio\cr
Type: \tab Package\cr
Version: \tab 1.2.4\cr
Date: \tab 2025-05-13\cr
License: \tab GPL 3.0\cr
LazyLoad: \tab yes\cr
}
The datasets are: Nuclei, Groups, CondA, CondB, Scores, and Dyad

The main functions are: boxWhisker, completeData, F_select, intPlot, ldaPlot, loadings, meanCent, percentMax, permuteLDA, power, ppca_mdb, zTrans, binomPower, h2Estimate, and plotBinomPower.

Type ?<object> to learn more about these objects, e.g. ?Nuclei

Type ?<function> to see examples of the function's use, e.g. ?FSelect
}

\author{
	Samuel V Scarpino
	Maintainer: Samuel V Scarpino <scarpino@utexas.edu>
}

\references{
Collyer M, Adams D. (2007) Analysis of Two - State Multivariate Phenotypic Change in Ecological Studies. Ecology: 88(3) 683 - 692.

Costanza M, Afifi A. (1979) Comparison of Stopping Rules in Forward Stepwise Discriminant Analysis. Journal of the American Statistical Association: pp. 777 - 78

Crews D, Gillette R, Scarpino SV, Manikkam M, Savenkova MI, Skinner MK. (2012)
Epigenetic Transgenerational Alterations to Stress Response in Brain Gene Networks and 
Behavior. Proc. Natl. Acad. Sci. USA: 109(23) 9143 - 9148.

Davies SW, Scarpino SV, Pongwarin T, Scott J, Matz MV. (2015) Estimating Trait Heritability in Highly Fecund Species. G3: Genes| Genomes| Genetics: 5(12) 2639 - 45.

Habbema J, Hermans J. (1977) Selection of Variables in Discriminant Analysis by F-Statistics and Error Rate. Technometrics: 19(4) 487 - 493.

Jennrich R. (1977) Stepwise discriminant analysis, volume 3. New York Wiley Sons.

Roweis S. (1997) EM algorithms for PCA and sensible PCA. Neural Inf. Proc. Syst.: 10
626 - 632.

Stacklies W, Redestig H, Scholz M, Walther D, Selbig J. (2007) pcaMethods - a Bioconductor package providing PCA methods for incomplete data. Bioinformatics: 23 1164 - 1167.

Troyanskaya O, Cantor M, Sherlock G, Brown P, Hastie T, Tibshirani R, Botstein D, Altman R. (2001) Missing value estimation methods for DNA microarrays. Bioinformatics: 17(6) 520 - 5252.
}

\seealso{

\code{\link[pcaMethods:pcaMethods]{pcaMethods}}
}
