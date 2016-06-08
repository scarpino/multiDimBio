#' Function to impute missing data.
#' 
#' This function imputes missing data using a probabilistic principle component
#' analysis framework and is a wrapper around functions implemented in the
#' pcaMethods package (Stacklies et al. 2007), was proposed by Troyanskaya et
#' al 2001 and is based on methods developed in Roweis 1997.
#' 
#' 
#' @param data a (non-empty) numeric matrix of data values.
#' @param n_pcs a (non-empty) numeric value indicating the desired number of
#' principle component axes.
#' @param cut.trait a number indicating the maximum proportion of missing
#' traits before an individual is removed from data.  A value of 1 will not
#' remove any individuals and 0 will remove them all.
#' @param cut.ind a number indicating the maximum proportion of individuals
#' missing a trait score before that trait is removed from data.  A value of 1
#' will not remove any traits and 0 will remove them all.
#' @param show.test a logical statement indicating whether a diagnostic plot of
#' the data imputation should be returned.
#' @return Returns a list with two entries.  \item{complete_dat}{ an object of
#' class matrix with missing values imputed using a probabilistic principle
#' component framework. } \item{plots}{ a list of plots stored as grid plots. }
#' @seealso \code{\link[pcaMethods:pcaMethods]{pcaMethods}}, \code{\link{pca}}
#' @references Roweis S (1997). EM algorithms for PCA and sensible PCA. Neural
#' Inf. Proc. Syst., 10, 626 - 632.
#' 
#' Stacklies W, Redestig H, Scholz M, Walther D, Selbig J (2007). pcaMethods -
#' a Bioconductor package providing PCA methods for incomplete data.
#' Bioinformatics, 23, 1164 - 1167.
#' 
#' Troyanskaya O, Cantor M, Sherlock G, Brown P, Hastie T, Tibshirani R,
#' Botstein D, Altman R (2001). Missing value estimation methods for DNA
#' microarrays. Bioinformatics, 17(6), 520 - 5252.
#' @examples
#' 
#' data(Nuclei)
#' npcs<-floor(ncol(Nuclei)/5)
#' 
#' length(which(is.na(Nuclei))==TRUE)
#' 
#' dat.comp<-completeData(data = Nuclei, n_pcs = npcs)
#' 
#' length(which(is.na(dat.comp))==TRUE)
#' 
#' @export completeData
completeData <- function(data, n_pcs, cut.trait = 0.5, cut.ind = 0.5, 
    show.test = TRUE) {
    
    # removing traits
    cut.trait <- round(nrow(data) * cut.trait)
    na.mat.tr <- is.na(data)
    na.sum.tr <- colSums(na.mat.tr)
    rm.na.tr <- which(na.sum.tr > cut.trait)
    
    if (length(rm.na.tr) > 0) {
        data <- data[ ,-rm.na.tr]
        cat("Traits removed: ", rm.na.tr, "\n", "\n")
    }
    
    # removing individuals
    cut.ind <- round(ncol(data) * cut.ind)
    na.mat.in <- is.na(data)
    na.sum.in <- rowSums(na.mat.in)
    rm.na.in <- which(na.sum.in > cut.ind)
    
    if (length(rm.na.in) > 0) {
        data <- data[-rm.na.in, ]
        cat("Individuals removed:", rm.na.in, "\n", "\n")
    }
    
    # ppca
    ppca <- pca(data, n_pcs = n_pcs, method = "ppca", center = TRUE, 
        scale = "vector")
    
    # Impute missing data
    imp1 <- completeObs(ppca)
    plots_ret <- list()
    
    if (show.test == TRUE) 
        {
            na.mat <- is.na(data)
            na.sum <- rowSums(na.mat)
            rm.NA <- which(na.sum > 0)
            
            if (length(rm.NA) > 0) 
                {
                  data.complete <- as.matrix(data[-rm.NA, ])
                }  #end if length(rm.NA)
            
            missing <- sum(na.sum)/(ncol(data) * nrow(data))
            censor <- round((ncol(data.complete) * nrow(data.complete)) * 
                missing)
            dataC <- matrix(data.complete, ncol = 1)
            make.NA <- sample(1:nrow(dataC), censor, replace = FALSE)
            dataC[make.NA, ] <- NA
            data.c.NA <- matrix(dataC, ncol = ncol(data.complete))
            
            # ppca
            ppca.c.NA <- pca(data.c.NA, n_pcs = n_pcs, method = "ppca", 
                center = TRUE, scale = "vector")
            
            # Impute missing data
            imp2 <- completeObs(ppca.c.NA)
           
            for (p in 1:ncol(imp2)) {
                layout(matrix(c(3, 1, 1, 3, 3, 1, 1, 3, 3, 2, 
                  2, 3, 3, 2, 2, 3), ncol = 4, byrow = TRUE))
                cols <- c("black", "red")
                ids <- rep(1, nrow(imp2))
                ids[which(is.na(data.c.NA[, p]) == TRUE)] <- 2
                plot(imp2[, p], data.complete[, p], main = colnames(data.complete)[p], 
                  xlab = "Imputed data", ylab = "Real data", pch = 16, 
                  col = cols[ids])
                hist((imp2[, p] - data.complete[, p])/mean(data.complete[, 
                  p]), breaks = 10, xlab = "Relative Error", main = colnames(data.complete)[p], 
                  col = "grey")
				plots_ret[[p]] <- recordPlot()
            }  #end for p
        }  #end if show.test==TRUE
    
    return(list("complete_dat" = imp1, "plots" = plots_ret))
    
}  #end function
