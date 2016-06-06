#' A function to visualize trait loadings onto discriminant function and
#' principle component axes
#' 
#' This function produces barplots representative of the contribution of a
#' particular trait or variable to either a discriminant function or principle
#' component axis.
#' 
#' 
#' @param DATA A (non-empty) numeric matrix with trait values
#' @param GROUPS A (non-empty)factor vector indicating the group membership of
#' each row in DATA
#' @param method An optional list indicating whether the results for a
#' principle component analysis, 'PCA', or linear discriminant analysis, 'LDA'
#' should be performed.
#' @return Outputs a .pdf file for each test listed in method.
#' @seealso \code{\link{pca}}, \code{\link{lda}}
#' @examples
#' 
#' data(Nuclei) 
#' data(Groups) 
#' Loadings(Nuclei, Groups, method=c("PCA", "LDA"))
#' 
#' @export Loadings
Loadings <- function(DATA, GROUPS, method = c("PCA", "LDA")) {
    
    results <- list()
    plots_ret <- list()
    
    if (sum(method == "PCA") > 0) 
        {
            nPCS = floor(ncol(DATA)/5)
            
            PPCA <- pca(DATA, nPcs = nPCS, method = "ppca", center = TRUE, 
                scale = "vector")
            
            OUT <- PPCA@loadings
            rownames(OUT) <- 1:nrow(OUT)
            NAMES <- data.frame(1:nrow(OUT), rownames(PPCA@loadings))
            colnames(NAMES) <- c("Number", "Trait")

            
            results[["Number_Trait_PCA"]] <- NAMES
            
            results[["Loadings"]] <- PPCA@loadings
            
            for (i in 1:nPCS) {
               
                title <- paste("PC", i, sep = "")
                barplot(abs(OUT[, i]), main = paste(title, "- Variance Explained = ", 
                  round(PPCA@R2[i], 3)), cex.names = 0.5)
                
                plots_ret[[paste(i, "PCA-Loadings.pdf", sep = "_")]] <- recordPlot()
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
     
      
            results[["Number_Trait_PCA"]] <- NAMES
            
            results[["Loadings"]] <- LDA$scaling
            
            for (j in 1:ncol(OUT)) {
                title <- paste("LD", j, sep = "")
                barplot(abs(OUT[, j]), main = title, cex.names = 0.5)
                plots_ret[[paste(j, "LDA-Loadings.pdf", sep = "_")]] <- recordPlot()
            }  #end for j
        }  #end if LDA
    return(list("results" = results, "plots" = plots_ret))
}  #end FUNCTION

