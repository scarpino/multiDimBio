#' A function to convert data into a z-score
#' 
#' This function converts the columns in a data matrix into z-scores.  The
#' score is computed by subracting each observation in a column from the column
#' mean and divding by the column standard deviation.  Each column is converted
#' independently of the others missing values are ignored in the calculation.
#' 
#' 
#' @param DATA A (non-empty) matrix with data values.  Columns should be
#' different traits and rows unique observations of those traits
#' @return Returns a matrix with the same dimensions as DATA.
#' @seealso \code{\link{PercentMax}}, \code{\link{MeanCent}}
#' @examples
#' 
#' data(Nuclei)
#' 
#' colMeans(Nuclei, na.rm=TRUE)
#' 
#' Nuclei.ZT<-ZTrans(Nuclei)
#' 
#' colMeans(Nuclei.ZT, na.rm=TRUE)
#' 
#' @export ZTrans
ZTrans <- function(DATA) {
    cmean <- apply(DATA, 2, mean, na.rm = TRUE)
    cmat <- matrix(cmean, ncol = ncol(DATA), nrow = nrow(DATA), 
        byrow = TRUE)
    muDAT <- DATA - cmat
    
    csd <- apply(DATA, 2, sd, na.rm = TRUE)
    cdmat <- matrix(csd, ncol = ncol(DATA), nrow = nrow(DATA), 
        byrow = TRUE)
    
    Z.SCORES <- muDAT/cdmat
    
    return(Z.SCORES)
}  #end FUNCTION

