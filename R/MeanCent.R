#' A function to scale data to mean 0
#' 
#' This function rescales the columns in a data matrix to have mean 0.  The
#' variance is not scaled and missing values are ignored in the calculation.
#' 
#' 
#' @param DATA A (non-empty) matrix with data values.  Columns should be
#' different traits and rows unique observations of those traits
#' @return Returns a matrix with the same dimensions as DATA.
#' @seealso \code{\link{ZTrans}}, \code{\link{PercentMax}}
#' @examples
#' 
#' data(Nuclei)
#' 
#' colMeans(Nuclei, na.rm=TRUE)
#' 
#' Nuclei.MC<-MeanCent(Nuclei)
#' 
#' colMeans(Nuclei.MC, na.rm=TRUE)
#' 
#' @export MeanCent
MeanCent <- function(DATA) {
    cmean <- apply(DATA, 2, mean, na.rm = TRUE)
    cmat <- matrix(cmean, ncol = ncol(DATA), nrow = nrow(DATA), 
        byrow = TRUE)
    muDAT <- DATA - cmat
    
    return(muDAT)
}  #end FUNCTION

