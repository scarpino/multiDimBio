MeanCent <- function(DATA) {
    cmean <- apply(DATA, 2, mean, na.rm = TRUE)
    cmat <- matrix(cmean, ncol = ncol(DATA), nrow = nrow(DATA), 
        byrow = TRUE)
    muDAT <- DATA - cmat
    
    return(muDAT)
}  #end FUNCTION

