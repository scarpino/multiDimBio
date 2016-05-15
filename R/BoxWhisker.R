#' A function to create a box and whisker plot by group ID
#' 
#' A function to create a box and whisker plot by group ID.
#' 
#' 
#' @param DATA a (non-empty) matrix of data values
#' @param GROUPS a (non-empty) vector of Group IDs with length equal to the
#' number of rows in DATA
#' @param palette A color palette for plotting.  The default is 'Paired.'  See
#' colorbrewer2.org for alternatives.
#' @return Saves a box-whisker plot of the data by group ID as a .pdf in the
#' working directory.
#' @examples
#' 
#' data(Nuclei)
#' data(Groups)
#' BoxWhisker(Nuclei, Groups)
#' 
#' #changing the color palette
#' 
#' BoxWhisker(Nuclei, Groups, palette='Set1')
#' 
#' @export BoxWhisker
BoxWhisker <- function(DATA, GROUPS, palette = "Paired") {
    
    message <- paste("Using", palette, sep = " ")
    
    cat(paste(message, "palette from colorbrewer2.org", sep = " "), 
        "\n", "\n")
    
    DAT <- as.matrix(DATA)
    SCORE <- matrix(DAT, ncol = 1)
    TRAIT <- colnames(DATA)
    TRAIT <- rep(TRAIT, each = nrow(DAT))
    GROUP <- rep(GROUPS, ncol(DAT))
    
    DAT.df <- data.frame(SCORE, TRAIT, GROUP)
    DAT.df$GROUP <- factor(DAT.df$GROUP)
    DAT.TRAIT <- factor(DAT.df$TRAIT)
    
    p <- ggplot(DAT.df, aes(TRAIT, SCORE))
    p.save <- p + geom_boxplot(aes(fill = GROUP), outlier.shape = NA) + 
        scale_fill_brewer(palette = palette) + theme(legend.position = "right", 
        legend.background = element_rect(fill = "#ffffffaa", colour = "black"), 
        panel.background = element_rect(fill = "white", colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 15), 
        axis.text.x = element_text(colour = "black", size = 8), 
        axis.title = element_text(colour = "black", size = 15), 
        legend.key = element_rect(fill = "white "), panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())
    
    timestamp <- as.character(as.integer(Sys.time()))
    ggsave(filename = paste(timestamp, "BoxWhisker.pdf", sep = "_"), 
        p.save)
}  #end FUNCTION

