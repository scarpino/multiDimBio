IntPlot <- function(Scores, Cov.A, Cov.B, pvalues = rep(1, 8), 
    int.pvalues = rep(1, 4)) {
    
    cat("It is best to use this function with PC scores", "\n", 
        "\n")
    
    if (ncol(Scores) > 4) {
        cat("Only the first 4 Score axes are plotted", "\n", "\n")
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
    
    pdf("Pvalue-Colors-InteractionPlots.pdf")
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
    dev.off()
    
    timestamp <- as.character(as.integer(Sys.time()))
    pdf(paste(timestamp, "Interaction Plots-LEGEND.pdf", sep = "-"))
    layout(matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE))
    
    for (i in 1:4) {
        interaction.plot(Cov.A, Cov.B, Scores[, i], legend = TRUE, 
            ylab = paste("Axis", i, "Score", sep = " "), xlab = "", 
            col = c(COL[i, 1], COL[i, 2]), lty = LTY[i], lwd = 2)
    }
    dev.off()
    
    pdf(paste(timestamp, "Interaction Plots.pdf", sep = "-"))
    layout(matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE))
    
    for (i in 1:4) {
        interaction.plot(Cov.A, Cov.B, Scores[, i], legend = FALSE, 
            ylab = paste("Axis", i, "Score", sep = " "), xlab = "", 
            col = c(COL[i, 1], COL[i, 2]), lty = LTY[i], lwd = 2)
    }
    dev.off()
    
}  #end FUNCTION

