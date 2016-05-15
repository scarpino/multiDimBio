partialF <- function(m.lda, GROUP, T_pm1) {
    GRS <- unique(GROUP)
    n1 <- length(which(GROUP == GRS[1]))
    n2 <- length(which(GROUP == GRS[2]))
    v <- n1 + n2 - 2
    p <- ncol(m.lda$means)
    T_p <- sum(m.lda$svd)
    T_pm1 <- T_pm1
    
    Fp <- (v - p + 1) * ((T_p - T_pm1)/(v + T_pm1))
    return(Fp)
}  #end partialF

