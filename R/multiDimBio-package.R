

#' Treatment condition for animals contained in the data set Nuclei
#' 
#' Animals measured in the Nuclei data set were either from linneages exposed
#' to the fungicide Vinclozolin (Vinclozolin) or not (Control).
#' 
#' 
#' @name CondA
#' @docType data
#' @format A factor vector indicating which treatment group the individuals in
#' Nuclei belong to.
#' @references Crews, D, R Gillette, SV Scarpino, M Manikkam, MI Savenkova, MK
#' Skinner. 2012. Epigenetic Transgenerational Alterations to Stress Response
#' in Brain Gene Networks and Behavior. Proc. Natl. Acad. Sci. USA. 109 (23).
#' 9143 - 9148.
#' @source The data are provided courtesy of David Crews at the University of
#' Texas at Austin.
#' @examples
#' 
#' data(CondA)
#' 
NULL





#' Stress condition for animals contained in the data set Nuclei
#' 
#' Animals measured in the Nuclei data set were either subjected to chronic
#' restraint stress (stress) or not (control).
#' 
#' 
#' @name CondB
#' @docType data
#' @format A factor vector indicating which stress group the individuals in
#' Nuclei belong to.
#' @references Crews, D, R Gillette, SV Scarpino, M Manikkam, MI Savenkova, MK
#' Skinner. 2012. Epigenetic Transgenerational Alterations to Stress Response
#' in Brain Gene Networks and Behavior. Proc. Natl. Acad. Sci. USA. 109 (23).
#' 9143 - 9148.
#' @source The data are provided courtesy of David Crews at the University of
#' Texas at Austin.
#' @examples
#' 
#' data(CondB)
#' 
NULL





#' Housing dyad for animals contained in the data set Nuclei
#' 
#' Animals measured in the Nuclei data set were housed in dyads with one
#' individual from the Vinclozolin line and one from the control line housed
#' together.  Each dyad was either stressed or not stressed.
#' 
#' 
#' @name Dyad
#' @docType data
#' @format A factor vector indicating which housing dyad the individuals in
#' Nuclei are in.
#' @references Crews, D, R Gillette, SV Scarpino, M Manikkam, MI Savenkova, MK
#' Skinner. 2012. Epigenetic Transgenerational Alterations to Stress Response
#' in Brain Gene Networks and Behavior. Proc. Natl. Acad. Sci. USA. 109 (23).
#' 9143 - 9148.
#' @source The data are provided courtesy of David Crews at the University of
#' Texas at Austin.
#' @examples
#' 
#' data(Dyad)
#' 
NULL





#' The group ID for animals contained in the data set Nuclei
#' 
#' Animals measured in the Nuclei data set belong to one of four groups
#' determined by their linneage (Vinclozolin or Control) and their stress
#' treatment (Stressed or Non-Stressed).
#' 
#' 
#' @name Groups
#' @docType data
#' @format A factor vector indicating which group the individuals in Nuclei are
#' in.
#' @references Crews, D, R Gillette, SV Scarpino, M Manikkam, MI Savenkova, MK
#' Skinner. 2012. Epigenetic Transgenerational Alterations to Stress Response
#' in Brain Gene Networks and Behavior. Proc. Natl. Acad. Sci. USA. 109 (23).
#' 9143 - 9148.
#' @source The data are provided courtesy of David Crews at the University of
#' Texas at Austin.
#' @examples
#' 
#' data(Groups)
#' 
NULL





#' A Package for the Design, Analysis, and Visualization of Systems Biology
#' Experiments
#' 
#' multiDimBio is a package designed to support a systems biology research
#' program from inception through publication.  It focuses on dimension
#' reduction approaches to detect patterns in complex, multivariate
#' experimental data and places an emphasis on informative visualizations.  The
#' goal for this project is to create a package that will evolve over time,
#' thereby remaining relevant and reflective of current methods and techniques.
#' As a result, we encourage suggested additions to the package, both
#' methodological and graphical.
#' 
#' \tabular{ll}{ Package: \tab multiDimBio\cr Type: \tab Package\cr Version:
#' \tab 1.2.3\cr Date: \tab 2025-05-13\cr License: \tab GPL 3.0\cr LazyLoad:
#' \tab yes\cr } The datasets are: Nuclei, Groups, CondA, CondB, Scores, and
#' Dyad
#' 
#' The main functions are: boxWhisker, completeData, F_select, intPlot,
#' ldaPlot, loadings, meanCent, percentMax, permuteLDA, power, ppca_mdb,
#' zTrans, binomPower, h2Estimate, and plotBinomPower.
#' 
#' Type ?<object> to learn more about these objects, e.g. ?Nuclei
#' 
#' Type ?<function> to see examples of the function's use, e.g. ?FSelect
#' 
#' @name multiDimBio-package
#' @aliases multiDimBio-package multiDimBio
#' @docType package
#' @author Samuel V Scarpino Maintainer: Samuel V Scarpino
#' <scarpino@@utexas.edu>
#' @seealso
#' 
#' \code{\link[pcaMethods:pcaMethods]{pcaMethods}}
#' @references Collyer M, Adams D. (2007) Analysis of Two - State Multivariate
#' Phenotypic Change in Ecological Studies. Ecology: 88(3) 683 - 692.
#' 
#' Costanza M, Afifi A. (1979) Comparison of Stopping Rules in Forward Stepwise
#' Discriminant Analysis. Journal of the American Statistical Association: pp.
#' 777 - 78
#' 
#' Crews D, Gillette R, Scarpino SV, Manikkam M, Savenkova MI, Skinner MK.
#' (2012) Epigenetic Transgenerational Alterations to Stress Response in Brain
#' Gene Networks and Behavior. Proc. Natl. Acad. Sci. USA: 109(23) 9143 - 9148.
#' 
#' Davies SW, Scarpino SV, Pongwarin T, Scott J, Matz MV. (2015) Estimating
#' Trait Heritability in Highly Fecund Species. G3: Genes| Genomes| Genetics:
#' 5(12) 2639 - 45.
#' 
#' Habbema J, Hermans J. (1977) Selection of Variables in Discriminant Analysis
#' by F-Statistics and Error Rate. Technometrics: 19(4) 487 - 493.
#' 
#' Jennrich R. (1977) Stepwise discriminant analysis, volume 3. New York Wiley
#' Sons.
#' 
#' Roweis S. (1997) EM algorithms for PCA and sensible PCA. Neural Inf. Proc.
#' Syst.: 10 626 - 632.
#' 
#' Stacklies W, Redestig H, Scholz M, Walther D, Selbig J. (2007) pcaMethods -
#' a Bioconductor package providing PCA methods for incomplete data.
#' Bioinformatics: 23 1164 - 1167.
#' 
#' Troyanskaya O, Cantor M, Sherlock G, Brown P, Hastie T, Tibshirani R,
#' Botstein D, Altman R. (2001) Missing value estimation methods for DNA
#' microarrays. Bioinformatics: 17(6) 520 - 5252.
NULL





#' Brain activity in 14 brain regions for 71 individuals
#' 
#' The activity in 14 brain nuclei were measured in rats that were in one of
#' four groups: 1) Non-stressed, Control 2) Stressed, Control 3) Non-stressed,
#' Vinclozolin 4) Stressed, Vinclozolin
#' 
#' Two different cohorts of male rats of the F3 generation of Vinclozolin
#' (Vinclozolin-Lineage) and Vehicle Control (Control-Lineage) Lineages
#' produced at Washington State University are shipped to the University of
#' Texas on the day after weaning. Rats are randomly pair-housed (one
#' Control-Lineage and one Vinclozolin-Lineage animal) and remain in these
#' dyads throughout the duration of the study. Half of the dyads are randomly
#' chosen to receive chronic restraint stress (CRS) treatment for 6 hours daily
#' for 21 consecutive days commencing 1 hr after lights off.  Activity in 14
#' brain nuclei were measured at the end of the study.
#' 
#' @name Nuclei
#' @docType data
#' @format A numeric matrix with 71 individuals as rows and the activity of 14
#' brain nuclei as columns.  NAs indicate missing data.
#' @references Crews, D, R Gillette, SV Scarpino, M Manikkam, MI Savenkova, MK
#' Skinner. 2012. Epigenetic Transgenerational Alterations to Stress Response
#' in Brain Gene Networks and Behavior. Proc. Natl. Acad. Sci. USA. 109 (23).
#' 9143 - 9148.
#' @source The data are provided courtesy of David Crews at the University of
#' Texas at Austin.
#' @examples
#' 
#' data(Nuclei)
#' 
NULL





#' Principle component scores based on the data in Nuclei
#' 
#' Principle component scores were computed using PPCA for the data set Nuclei.
#' 
#' 
#' @name Scores
#' @docType data
#' @format A numeric matrix with 4 columns and the same number of rows as
#' Nuclei.  There are no missing values.
#' @references Crews, D, R Gillette, SV Scarpino, M Manikkam, MI Savenkova, MK
#' Skinner. 2012. Epigenetic Transgenerational Alterations to Stress Response
#' in Brain Gene Networks and Behavior. Proc. Natl. Acad. Sci. USA. 109 (23).
#' 9143 - 9148.
#' @source The data are provided courtesy of David Crews at the University of
#' Texas at Austin.
#' @examples
#' 
#' data(Scores)
#' 
#' data(Nuclei)
#' 
#' SCORES<-PPCA(Nuclei)@scores
#' 
NULL



