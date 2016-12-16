#' A function to plot the results of a binomPower run
#' 
#' A function to plot the results of a binomPower run.
#' 
#' 
#' @param datPlotBig a (non-empty) matrix of data values, with columns trueTau,
#' ndads, trueTau2
#' @param params a (non-empty) matrix of parameter values, with columns mm and
#' vv.
#' @return Returns a list of two plots of the binomPower analysis results.
#' @examples
#' 
#' #not run
#' 
#' @export plotBinomPower
plotBinomPower <- function(datPlotBig, params) {
    
    if (length(unique(params$mm)) > 1) {
        stop("function does not work with more than one mean offspring number")
    }
    
    if (length(unique(params$vv)) > 1) {
        stop("function does not work with more than one var in offspring number")
    }
    
    # solving ggplot2 binding error
    ndads <- NULL
      
    # preping data setss
    datPlotBig$trueTau <- as.factor(datPlotBig$trueTau)
    datPlotBig$ndads <- round(as.numeric(as.character(datPlotBig$ndads)))
    datPlotBig$ndads <- round(datPlotBig$ndads, 1)
    datPlotBig$ndads <- as.factor(datPlotBig$ndads)
    
    # colors
    cols <- colorRampPalette(c("#FFF7FB", "#74A9CF", "#023858"), 
        space = "Lab")
    pal <- cols(length(unique(datPlotBig$ndads)))
    
    ###### plotting#
    
    plots_ret <- list()
    
    # values h2 > 0
    use <- which(datPlotBig$trueTau2 > 0)
    dat.plot <- datPlotBig[use, ]
    
    # transforming variances to heritability
    h2 <- 4 * (dat.plot$trueTau2/(dat.plot$trueTau2 + (pi/sqrt(3))^2))
    auc <- dat.plot$trueP/(dat.plot$trueP + dat.plot$falseP)
    dat.plot <- data.frame(dat.plot, h2, auc)
    
    p <- ggplot(dat.plot, aes(h2, auc, colour = ndads, group = ndads))
    p.save <- p + geom_line(size = 1.5) + scale_colour_manual(values = pal) + 
        xlab("h2") + ylab("auc") + theme(legend.position = "right", 
        legend.key = element_rect(fill = "gray"), legend.background = element_rect(fill = "#ffffffaa", 
            colour = "black"), panel.background = element_rect(fill = "gray", 
            colour = "black"), axis.text.y = element_text(colour = "black", 
            size = 15), axis.text.x = element_text(colour = "black", 
            size = 20), axis.title = element_text(colour = "black", 
            size = 20), 
        panel.grid.minor = element_line(colour = "#00000050", 
            linetype = 3), panel.grid.major = element_line(colour = "#00000060", 
            linetype = 3)) + scale_y_continuous(expand = c(0.005, 
        0.005)) + scale_x_continuous(expand = c(0, 0.003), limits = c(0, 
        max(dat.plot$h2)))
    
    plots_ret[["powerCurveBinom"]] <- p.save
    
    # values h2 == 0
    use <- which(datPlotBig$trueTau2 == 0)
    dat.plot <- datPlotBig[use, ]
    auc <- dat.plot$trueP/(dat.plot$trueP + dat.plot$falseP)
    rownames(dat.plot) <- NULL
    dat.plot <- data.frame(dat.plot, h2, auc)
    vals <- by(dat.plot[, "auc"], dat.plot[, "ndads"], mean, na.rm = TRUE)
    dat.plot <- data.frame(as.factor(as.numeric(names(vals))), 
        unlist(vals)[1:length(vals)])
    colnames(dat.plot) <- c("ndads", "auc")
    
    p1 <- ggplot(dat.plot, aes(x = ndads, y = auc, fill = ndads))
    p1.save <- p1 + geom_bar(stat = "identity") + scale_fill_manual(values = pal) + 
        xlab("number of dads") + ylab("auc (true h2 = 0)") + theme(legend.position = "right", 
        legend.key = element_rect(fill = "gray"), legend.background = element_rect(fill = "#ffffffaa", 
            colour = "black"), panel.background = element_rect(fill = "gray", 
            colour = "black"), axis.text.y = element_text(colour = "black", 
            size = 15), axis.text.x = element_text(colour = "black", 
            size = 20), axis.title = element_text(colour = "black", 
            size = 20), 
        panel.grid.minor = element_line(colour = "#00000050", 
            linetype = 3), panel.grid.major = element_line(colour = "#00000060", 
            linetype = 3)) + scale_y_continuous(expand = c(0, 
        0), limits = c(0, 1)) + scale_x_discrete(expand = c(0.01, 
        0.01))
    
    plots_ret[["powerTrue_h2_equal_0_Binom"]] <- p1.save
    
    return(plots_ret)
}
