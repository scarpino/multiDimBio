#' A function to visualize the results of a discriminant analysis
#' 
#' The function takes as input the traits and group IDs and will perform a
#' discriminate function analysis and visualize the results. For the pair-wise
#' comparison of groups we use density histograms with points along the x-axis
#' denoting the actual data, Figure 3 For multi-group comparisons we plot a
#' bivariate scatter for all pairwise combinations of discriminate axes. The
#' color of plotting symbols can be altered using the palette argument and the
#' axes comparisons (with max n = number of groups - 1).
#' 
#' 
#' @param Data A (non-empty), numeric matrix of data values
#' @param Groups A (non-empty), vector indicating group membership.
#' Length(unique(Group))==2
#' @param palette A color palette for plotting.  The default is 'Paired.'  See
#' colorbrewer2.org for alternatives.
#' @param axes A numeric vector describing which axes to compare.  For example,
#' axes=c(1,2) will on produce a single plot comparing the first and second
#' axis.
#' @return Returns a list of ggplot2 plots.
#' @seealso \code{\link{lda}}
#' @examples
#' 
#' data(Nuclei)
#' data(Groups) 
#' ldaPlot(Nuclei, Groups, palette='BrBG', axes=c(1,2,2,3,1,3))
#' 
#' @export ldaPlot
ldaPlot <- function(Data, Groups, palette = "BrBG", axes = c(1, 
    2, 2, 3, 1, 3)) {
    LD1 <- NULL
    LD2 <- NULL
    
    if (min(Data, na.rm = TRUE) < 1e-04) {
        TOL <- min(Data, na.rm = TRUE)
    } else {
        TOL <- 1e-04
    }
    
    plots_ret <- list()
    
    Groups <- as.factor(Groups)
    Groups <- as.numeric(Groups)
    
    if (length(unique(Groups)) > 2) 
        {
            
            if (length(unique(Groups)) > 6) 
                {
                  stop("Will not work with more than 6 groups", 
                    "\n", "\n")
                }  #end if >6
            
            ld1 <- lda(Groups ~ Data, tol = TOL)
            pm1 <- predict(ld1)
            GC <- as.factor(pm1$class)
            pm1df <- data.frame(GC, pm1$x)

            
            AX <- matrix(axes, ncol = 2, byrow = TRUE)
            
            for (i in 1:nrow(AX)) {
                AX1 <- AX[i, 1] + 1
                AX2 <- AX[i, 2] + 1
                pm1df.i <- pm1df[, c(1, AX1, AX2)]
                colnames(pm1df.i) <- c("GC", "LD1", "LD2")
                # pi<-ggplot(pm1df,aes(pm1df[,AX1],pm1df[,AX2]))
                pi <- ggplot(pm1df.i, aes(LD1, LD2))
                pi.save <- pi + geom_point(aes(colour = GC, shape = GC), 
                  size = 3.5) + scale_colour_brewer(palette = palette) + 
                  scale_x_continuous(name = colnames(pm1df)[AX1]) + 
                  scale_y_continuous(name = colnames(pm1df)[AX2]) + 
                  theme(legend.position = "right", legend.background = element_rect(fill = "#ffffffaa", 
                    colour = "black"), panel.background = element_rect(fill = "white", 
                    colour = "black"), axis.text = element_text(colour = "black", 
                    size = 15), axis.title = element_text(colour = "black", 
                    size = 15), legend.key = element_rect(fill = "white "), 
                    panel.grid.minor = element_blank(), panel.grid.major = element_blank())
                
                
                plots_ret[["ldaPlot-AllGroups"]] <- pi.save

            }  #end for i
        }  #end if unique(groups)>2
    
    # 1 comparison matrix
    ngroups <- length(unique(Groups))
    comp <- makeCompMat(ngroups)
    
    for (j in 1:nrow(comp)) {
        use.j <- which(Groups == comp[j, 1] | Groups == comp[j, 
            2])
        
        ld.j <- lda(Groups[use.j] ~ Data[use.j, ], tol = TOL)
        pm.j <- predict(ld.j)
        GC.j <- as.factor(pm.j$class)
        pmdf.j <- data.frame(GC.j, pm.j$x)
        colnames(pmdf.j) <- c("GC", "LD1")
        
        p.j <- ggplot(pmdf.j, aes(x = LD1))
        
        psave.j <- p.j + geom_density(aes(fill = GC)) + scale_fill_manual(values = c("#8C510A66", 
            "#01665E66")) + geom_point(aes(y = 0, x = LD1, shape = GC), 
            size = 3.5) + geom_point(aes(y = 0, x = LD1, colour = GC, 
            shape = GC)) + scale_colour_manual(values = c("#8C510A", 
            "#01665E")) + labs(x = "Linear Discriminant 1 Score", 
            y = "Density") + theme(legend.position = "right", 
            legend.background = element_rect(fill = "#ffffffaa", 
                colour = "black"), panel.background = element_rect(fill = "white", 
                colour = "black"), axis.text = element_text(colour = "black", 
                size = 15), axis.title = element_text(colour = "black", 
                size = 15), legend.key = element_rect(fill = "white "), 
            panel.grid.minor = element_blank(), panel.grid.major = element_blank())
        
        timestamp <- as.character(as.integer(Sys.time()))
        plots_ret[[paste(comp[j, 1], comp[j, 2], "Group LDA Plot", sep = " ")]] <- psave.j
    }  #end for j
    return(plots_ret)
}  #end FUNCTION

