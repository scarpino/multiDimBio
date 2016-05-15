CompleteData <- function(DATA, NPCS, cut.trait = 0.5, cut.ind = 0.5, 
    show.test = TRUE) {
    
    # removing traits
    cut.trait <- round(nrow(DATA) * cut.trait)
    NA.mat.tr <- is.na(DATA)
    NA.sum.tr <- colSums(NA.mat.tr)
    rm.NA.tr <- which(NA.sum.tr > cut.trait)
    
    if (length(rm.NA.tr) > 0) {
        DATA <- DATA[, -rm.NA.tr]
        cat("Traits removed: ", rm.NA.tr, "\n", "\n")
    }
    
    # removing individuals
    cut.ind <- round(ncol(DATA) * cut.ind)
    NA.mat.in <- is.na(DATA)
    NA.sum.in <- rowSums(NA.mat.in)
    rm.NA.in <- which(NA.sum.in > cut.ind)
    
    if (length(rm.NA.in) > 0) {
        DATA <- DATA[-rm.NA.in, ]
        cat("Individuals removed:", rm.NA.in, "\n", "\n")
    }
    
    # PPCA
    PPCA <- pca(DATA, nPcs = NPCS, method = "ppca", center = TRUE, 
        scale = "vector")
    
    # Impute missing data
    imp1 <- completeObs(PPCA)
    
    if (show.test == TRUE) 
        {
            NA.mat <- is.na(DATA)
            NA.sum <- rowSums(NA.mat)
            rm.NA <- which(NA.sum > 0)
            
            if (length(rm.NA) > 0) 
                {
                  DATA.complete <- as.matrix(DATA[-rm.NA, ])
                }  #end if length(rm.NA)
            
            missing <- sum(NA.sum)/(ncol(DATA) * nrow(DATA))
            censor <- round((ncol(DATA.complete) * nrow(DATA.complete)) * 
                missing)
            DATAC <- matrix(DATA.complete, ncol = 1)
            make.NA <- sample(1:nrow(DATAC), censor, replace = FALSE)
            DATAC[make.NA, ] <- NA
            DATA.c.NA <- matrix(DATAC, ncol = ncol(DATA.complete))
            
            # PPCA
            PPCA.c.NA <- pca(DATA.c.NA, nPcs = NPCS, method = "ppca", 
                center = TRUE, scale = "vector")
            
            # Impute missing data
            imp2 <- completeObs(PPCA.c.NA)
            
            for (p in 1:ncol(imp2)) {
                pdf(paste(colnames(DATA.complete)[p], "ImputeTest.pdf", 
                  sep = "_"))
                layout(matrix(c(3, 1, 1, 3, 3, 1, 1, 3, 3, 2, 
                  2, 3, 3, 2, 2, 3), ncol = 4, byrow = TRUE))
                COLS <- c("black", "red")
                IDS <- rep(1, nrow(imp2))
                IDS[which(is.na(DATA.c.NA[, p]) == TRUE)] <- 2
                plot(imp2[, p], DATA.complete[, p], main = colnames(DATA.complete)[p], 
                  xlab = "Imputed Data", ylab = "Real Data", pch = 16, 
                  col = COLS[IDS])
                hist((imp2[, p] - DATA.complete[, p])/mean(DATA.complete[, 
                  p]), breaks = 10, xlab = "Relative Error", main = colnames(DATA.complete)[p], 
                  col = "grey")
                dev.off()
            }  #end for p
        }  #end if show.test==TRUE
    
    return(imp1)
    
}  #end FUNCTION

