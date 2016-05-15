#' A function to perform a probabilistic principle component analysis
#' 
#' Performs a probabilistic principle component analysis using the function
#' 'pca' in the package'pcaMethods'
#' 
#' In PPCA an Expectation Maximization (EM) algorithm is used to fit a Gaussian
#' latent variable model ( Tippping and Bishop (1999)). A latent variable model
#' seeks to relate an observed vector of data to a lower dimensional vector of
#' latent (or unobserved) variables, an approach similar to a factor analysis.
#' Our implementation is a wrapper around the pcaMethods functions ppca and
#' svdimpute (Stacklies et al. (2007)) and is included mainly for convience.
#' The method used in pca was adapted from Roweis (1997) and a Matlab script
#' developed by Jakob Verbeek.
#' 
#' @param Data A (non-empty), numeric matrix of data values
#' @param nPCs The number of resulting principle component axes.  nPCs must be
#' less than or equal to the number of columns in Data.
#' @param CENTER A logical statement indicating whether data should be centered
#' to mean 0, TRUE, or not, FALSE.
#' @param SCALE A character string indicating which method should be used to
#' scale the variances.  The default setting is 'vector.'
#' @return Returns an object of class 'pcaRes.'  See documentation in the
#' package code\link[pcaMethods:pcaMethods]{ pcaMethods}
#' @seealso \code{\link[pcaMethods:pcaMethods]{pcaMethods}}, \code{\link{pca}}
#' @references Roweis S (1997). EM algorithms for PCA and sensible PCA. Neural
#' Inf. Proc. Syst., 10, 626 - 632.
#' 
#' Stacklies W, Redestig H, Scholz M, Walther D, Selbig J (2007). pcaMethods -
#' a Bioconductor package providing PCA methods for incomplete data.
#' Bioinformatics, 23, 1164 - 1167.
#' 
#' Tippping M, Bishop C (1999). Probabilistic Principle Componenet Analysis.
#' Journal of the Royal Statistical Society. Series B (Statistical
#' Methodology), 61(3), 611 - 622.
#' @examples
#' 
#' 	data(Nuclei)
#' 	PPCA1<-PPCA(Nuclei, nPCs=2, CENTER=TRUE, SCALE='vector')
#' 	Scores1<-PPCA1@scores
#' 	
#' @export PPCA
PPCA <- function(Data, nPCs = 4, CENTER = TRUE, SCALE = "vector") {
    PC <- pca(Data, nPcs = nPCs, method = "ppca", center = CENTER, 
        scale = SCALE)
    return(PC)
}  #end Function
