PercentMax <- function(DATA) {
    cmax <- apply(DATA, 2, max, na.rm = TRUE)
    cmat <- matrix(cmax, ncol = ncol(DATA), nrow = nrow(DATA), 
        byrow = TRUE)
    PMX.DATA <- DATA/cmat
    return(PMX.DATA)
}  #end FUNCTION

