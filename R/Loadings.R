Loadings <- function(DATA, GROUPS, method = c("PCA", "LDA")) {
    
    if (sum(method == "PCA") > 0) 
        {
            nPCS = floor(ncol(DATA)/5)
            
            PPCA <- pca(DATA, nPcs = nPCS, method = "ppca", center = TRUE, 
                scale = "vector")
            
            OUT <- PPCA@loadings
            rownames(OUT) <- 1:nrow(OUT)
            NAMES <- data.frame(1:nrow(OUT), rownames(PPCA@loadings))
            colnames(NAMES) <- c("Number", "Trait")
            cat("Number to trait map, saved as a .csv in the working dir.", 
                "\n", "\n")
            print(NAMES)
            
            timestamp <- as.character(as.integer(Sys.time()))
            write.csv(NAMES, paste(timestamp, "Number_Trait_PCA.csv"), 
                row.names = FALSE)
            
            timestamp <- as.character(as.integer(Sys.time()))
            write.csv(PPCA@loadings, paste(timestamp, "PCA Loadings.csv"))
            
            timestamp <- as.character(as.integer(Sys.time()))
            for (i in 1:nPCS) {
                pdf(paste(timestamp, i, "PCA-Loadings.pdf", sep = "_"))
                title <- paste("PC", i, sep = "")
                barplot(abs(OUT[, i]), main = paste(title, "- Variance Explained = ", 
                  round(PPCA@R2[i], 3)), cex.names = 0.5)
                dev.off()
            }  #end for i
            
        }  #end if PCA
    
    if (sum(method == "LDA") > 0) 
        {
            
            if (min(DATA, na.rm = TRUE) < 1e-04) {
                TOL <- min(DATA, na.rm = TRUE)
            } else {
                TOL <- 1e-04
            }
            
            LDA <- lda(GROUPS ~ DATA, tol = TOL)
            
            OUT <- LDA$scaling
            rownames(OUT) <- 1:nrow(OUT)
            NAMES <- data.frame(1:nrow(OUT), rownames(LDA$scaling))
            colnames(NAMES) <- c("Number", "Trait")
            cat("Number to trait map, saved as a .csv in the working dir.", 
                "\n", "\n")
            print(NAMES)
            timestamp <- as.character(as.integer(Sys.time()))
            write.csv(NAMES, paste(timestamp, "Number_Trait_LDA.csv"), 
                row.names = FALSE)
            timestamp <- as.character(as.integer(Sys.time()))
            write.csv(LDA$scaling, paste(timestamp, "LDA Loadings.csv"))
            
            timestamp <- as.character(as.integer(Sys.time()))
            for (j in 1:ncol(OUT)) {
                pdf(paste(timestamp, j, "LDA-Loadings.pdf", sep = "_"))
                title <- paste("LD", j, sep = "")
                barplot(abs(OUT[, j]), main = title, cex.names = 0.5)
                dev.off()
            }  #end for j
        }  #end if LDA
    
}  #end FUNCTION

