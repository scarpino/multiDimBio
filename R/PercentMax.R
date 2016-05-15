#' A function to scale data to the percent of the maximum observed
#' 
#' This function rescales the columns in a data matrix to the percent of the
#' maximum observed value.  The variance is not scaled and missing values are
#' ignored in the calculation.
#' 
#' 
#' @param DATA A (non-empty) matrix with data values.  Columns should be
#' different traits and rows unique observations of those traits
#' @return Returns a matrix with the same dimensions as DATA.
#' @seealso \code{\link{ZTrans}}, \code{\link{MeanCent}}
#' @examples
#' 
#' data(Nuclei)
#' 
#' colMeans(Nuclei, na.rm=TRUE)
#' 
#' Nuclei.PM<-PercentMax(Nuclei)
#' 
#' colMeans(Nuclei.PM, na.rm=TRUE)
#' 
#' @export PercentMax
PercentMax <- function(DATA) {
    cmax <- apply(DATA, 2, max, na.rm = TRUE)
    cmat <- matrix(cmax, ncol = ncol(DATA), nrow = nrow(DATA), 
        byrow = TRUE)
    PMX.DATA <- DATA/cmat
    return(PMX.DATA)
}  #end FUNCTION

