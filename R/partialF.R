#' A function to compute partial F statistics
#' 
#' This is an internal function used in FSelect.  It can only be used for two
#' groups.  The partial F statistic is the additional contribution to the model
#' from adding one more trait.
#' 
#' 
#' @param m.lda An object of class 'lda'
#' @param GROUP A factor vector indicating group membership
#' @param T_pm1 The F statistic calculated for a discriminant analysis with
#' only the most informative trait.
#' @return Returns a partial F statistic
#' @seealso \code{\link{FSelect}}
#' @references Habbema J, Hermans J (1977). Selection of Variables in
#' Discriminant Analysis by F-Statistics and Error Rate. Technometrics, 19(4),
#' 487 - 493.
#' @examples
#' 
#' #Internal function used in FSelect
#' 
#' data(Nuclei)
#' data(Groups)
#' 
#' NPC<-floor(ncol(Nuclei)/5)
#' 
#' DAT.comp<-CompleteData(Nuclei, NPCS=NPC) 
#' Groups.use<-c(1,2) 
#' use.DAT<-which(Groups==Groups.use[1]|Groups==Groups.use[2])
#' 
#' DAT.use<-Nuclei[use.DAT,]
#' GR.use<-Groups[use.DAT]
#' 
#' traitA<-2
#' 
#' mlda<-MASS::lda(GR.use~DAT.use[,traitA])
#' 
#' F1<-partialF(mlda,GR.use,0)
#' 
#' traitB<-1
#' 
#' mlda2<-MASS::lda(GR.use~DAT.use[,c(traitA,traitB)])
#' 
#' partialF(mlda2,GR.use,F1)
#' 
#' 
#' @export partialF
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

