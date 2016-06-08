#' A function to create a box and whisker plot by group ID
#' 
#' A function to create a box and whisker plot by group ID.
#' 
#' 
#' @param data a (non-empty) matrix of data values
#' @param groups a (non-empty) vector of group IDs with length equal to the
#' number of rows in data
#' @param palette A color palette for plotting.  The default is 'Paired.'  See
#' colorbrewer2.org for alternatives.
#' @return Returns a box-whisker plot of the data by group ID.
#' @examples
#' 
#' data(Nuclei)
#' data(Groups)
#' boxWhisker(Nuclei, Groups)
#' 
#' #changing the color palette
#' 
#' boxWhisker(data = Nuclei, groups = Groups, palette = 'Set1')
#' 
#' @export boxWhisker
boxWhisker <- function(data, groups, palette = "Paired") {
    
    dat <- as.matrix(data)
    score <- matrix(dat, ncol = 1)
    trait <- colnames(data)
    trait <- rep(trait, each = nrow(dat))
    group <- rep(groups, ncol(dat))
    
    dat.df <- data.frame(score, trait, group)
    dat.df$group <- factor(dat.df$group)
    dat.trait <- factor(dat.df$trait)
    
    p <- ggplot(dat.df, aes(trait, score))
    p.save <- p + geom_boxplot(aes(fill = group), outlier.shape = NA) + 
        scale_fill_brewer(palette = palette) + theme(legend.position = "right", 
        legend.background = element_rect(fill = "#ffffffaa", colour = "black"), 
        panel.background = element_rect(fill = "white", colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 15), 
        axis.text.x = element_text(colour = "black", size = 8), 
        axis.title = element_text(colour = "black", size = 15), 
        legend.key = element_rect(fill = "white "), panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())
    
    return(p.save)
}  #end FUNCTION

