#' A function to visualize the results of a MANOVA
#' 
#' The function produces an interaction plot to demonstrate the results of a
#' MANOVA using the function interaction.plot.
#' 
#' 
#' @param Scores A (non-empty) numeric matrix of principle component scores or
#' raw data.
#' @param Cov.A A (non-empty) bivariate factor vector indicating the factor for
#' each row in Scores
#' @param Cov.B A (non-empty) bivariate factor vector indicating the factor for
#' each row in Scores
#' @param pvalues An optional vector of p values for each covariate across
#' Scores.  The length of pvalues must equal the number of columns in Scores
#' times 2.
#' @param int.pvalues An optional vector of p values for each interaction.  The
#' length of int.pvalues must equal the number of columns in Scores.
#' @return a list of plots stored as grid plots.
#' @seealso \code{\link{interaction.plot}}
#' @examples
#' 
#' data(Scores)
#' data(CondA)
#' data(CondB)
#' 
#' pvals<-c(0.03,0.6,0.05,0.07,0.9,0.2,0.5,0.3)
#' int.pvals<-c(0.3,0.45,0.5,0.12)
#' 
#' IntPlot(Scores,CondA,CondB,pvalues=pvals, int.pvalues=int.pvals)
#' 
#' @export IntPlot
IntPlot <- function(Scores, Cov.A, Cov.B, pvalues = rep(1, 8), 
    int.pvalues = rep(1, 4)) {
   
    
    if (ncol(Scores) > 4) {
        warning("Only the first 4 Score axes are plotted", "\n", "\n")
    }
    
    # line colors
    pvalues[pvalues > 0.05] <- 1
    pvalues[pvalues <= 0.05 & pvalues > 0.01] <- 2
    pvalues[pvalues <= 0.01 & pvalues > 0.005] <- 3
    pvalues[pvalues <= 0.005] <- 4
    colors <- c("gray", "#FDBB84", "#EF6548", "#990000")
    COLS <- colors[pvalues]
    COL <- matrix(COLS, ncol = 2, nrow = 4, byrow = TRUE)
    
    # line types for interaction
    int.pvalues[int.pvalues > 0.05] <- 5
    int.pvalues[int.pvalues <= 0.05] <- 1
    LTY <- int.pvalues
    
    plots_ret <- list()
    
    #interaction plot
    plot(rep(0, 4), 1:4, pch = 15, col = colors, cex = 5, xaxt = "n", 
        yaxt = "n", xlab = "", ylab = "", xlim = c(0, 0.5), ylim = c(0, 
            4), lwd = 2)
    segments(0, 0, 0.1, 0, lty = 3, lwd = 2)
    segments(0, 0.5, 0.1, 0.5, lty = 1, lwd = 2)
    text(0.1, 4, labels = "p<0.005", cex = 2)
    text(0.15, 3, labels = "0.01< p >0.005", cex = 2)
    text(0.15, 2, labels = "0.05< p >0.01", cex = 2)
    text(0.1, 1, labels = "p >0.05", cex = 2)
    text(0.27, 0.5, labels = "Interaction Significant", cex = 2)
    text(0.3, 0, labels = "Interaction Not Significant", cex = 2)
    
    plots_ret[["Pvalue-Colors-InteractionPlots"]] <- recordPlot()
    
   	#interaction plot legend
    layout(matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE))
    
    for (i in 1:4) {
        interaction.plot(Cov.A, Cov.B, Scores[, i], legend = TRUE, 
            ylab = paste("Axis", i, "Score", sep = " "), xlab = "", 
            col = c(COL[i, 1], COL[i, 2]), lty = LTY[i], lwd = 2)
    }
    plots_ret[["Interaction Plots-LEGEND"]] <- recordPlot()
    
    #interaction plots no colors
    layout(matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE))
    
    for (i in 1:4) {
        interaction.plot(Cov.A, Cov.B, Scores[, i], legend = FALSE, 
            ylab = paste("Axis", i, "Score", sep = " "), xlab = "", 
            col = c(COL[i, 1], COL[i, 2]), lty = LTY[i], lwd = 2)
    }
    plots_ret[["Interaction Plots"]] <- recordPlot()
    
    return(plots_ret)
}  #end function

