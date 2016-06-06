#' A function to determine whether two groups are in statistically different
#' locations in multivariate space See Collyer and Adams 2007
#' 
#' The function calculates the multivariate distance between two groups across
#' all traits and determines whether they differ signifcantly using a Monte
#' Carlo randomization test.  The Monte Carlo randomization creates a null
#' distribution by randomizing the residual deviation from the group mean
#' across all individuals.  This method controls for heteroscedasticity and was
#' designed by Collyer and Adams (2007) for use in analyzing data sets that
#' have sparse groups sizes relative to the number of traits.
#' 
#' Determining the statistical significance of a discriminate function analysis
#' along with performing that analysis on sparse data sets, e.g. many traits
#' observed on comparatively few individuals, is a challenge. Collyer and Adams
#' (2007) developed a Monte Carlo based algorithm for addressing both of those
#' issues. Briefly, the test uses the underlying Var/Cov structure of the data
#' and randomizes the group membership to calculate a null distribution. This
#' test simultaneously controls for heteroscedasticity, a common problem in
#' sparse data sets and allows the approximation of a p-value for the test. For
#' the original implementation and formulation of the method see Collyer and
#' Adams (2007) or http://www.public.iastate. edu/~dcadams/software.html.
#' Unlike the FSelect implementation, PermuteLDA will work properly with an
#' arbitrary number of groups. The time required to run the algorithm is
#' non-linear in the number of groups.
#' 
#' @param Data A (non-empty), numeric matrix of data values
#' @param Groups A (non-empty), vector indicating group membership.
#' @param NPerm The number of permutations used to generate the null
#' distribution.  The default is 100.
#' @param Missing.Data The method used to handle missing data.  The default,
#' 'Complete' will use completeData to impute missing data, setting
#' Missing.Data='Remove' will remove all individuals with missing data.
#' FSelect cannot handle missing data.
#' @return Returns a data frame with four columns and the number of groups
#' choose 2 rows.  Each row is a pairwise comparison between groups.  The
#' column 'Pr' is the p value to reject the null hypothesis of no difference (a
#' value in 'Pr' < 0.05 would result in rejecting the hypothesis that the two
#' groups are not different.  The column 'Distance' is the multivariate
#' distance between the two groups.
#' @seealso \code{\link{PermuteLDA}}
#' @references Collyer M, Adams D (2007). Analysis of Two - State Multivariate
#' Phenotypic Change in Ecological Studies. Ecology, 88(3), 683 - 692.
#' 
#' For an implementation of the original method coded in R see
#' http://www.public.iastate. edu/~dcadams/software.html.
#' @examples
#' 
#' 	data(Nuclei)
#' 	data(Groups)
#' 	PermuteLDA(Nuclei,Groups,50)
#' 
#' @export PermuteLDA
PermuteLDA <- function(Data, Groups, NPerm, Missing.Data = "Complete") {
    
    Groups <- as.factor(Groups)
    Groups <- as.numeric(Groups)
    
    # Missing data
    if (Missing.Data == "Complete" & length(which(is.na(Data) == 
        TRUE)) > 0) {
        warning("Missing data imputed using completeData to exlude missing data set Missing.Data=Remove", 
            "\n", "\n")
        Data <- completeData(Data, n_pcs = floor(ncol(Data)/5), 
            cut.trait = 1, cut.ind = 1)$complete_dat
    }
    
    if (Missing.Data == "Remove" & length(which(is.na(Data) == 
        TRUE)) > 0) {
        warning("Individuals with missing data removed, to impute missing data set Missing.Data=Complete", 
            "\n", "\n")
        rm <- na.omit(Data)
        Data <- Data[-attr(rm, "na.action"), ]
        Groups <- Groups[-attr(rm, "na.action")]
    }
    
    # 1 comparison matrix
    ngroups <- length(unique(Groups))
    comp <- makeCompMat(ngroups)
    
    
    # 2 comparisons
    
    results.i <- c()
    dists.i <- c()
    for (i in 1:nrow(comp)) {
        # selecting groups
        use.i <- which(Groups == comp[i, 1] | Groups == comp[i, 
            2])
        dat.i <- as.matrix(Data[use.i, ])
        group.i <- rep(1, nrow(dat.i))
        group2 <- which(Groups[use.i] == comp[i, 2])
        group.i[group2] <- -1
        
        # estimating parameters
        m.i <- lm(dat.i ~ group.i)
        
        # least-square mean
        obs.i <- by(m.i$fitted.value, group.i, colMeans)
        obs.i <- rbind(obs.i$"1", obs.i$"-1")
        
        # difference between groups
        obs.diff.i <- obs.i[1, ] - obs.i[2, ]
        dist.obs <- sqrt(t(obs.diff.i) %*% obs.diff.i)
        
        # expected values
        exp.i <- matrix(ncol = ncol(dat.i), nrow = nrow(dat.i))
        for (k in 1:ncol(exp.i)) {
            exp.i[which(group.i == 1), k] <- obs.i[1, k]
            exp.i[-which(group.i == 1), k] <- obs.i[2, k]
        }  #end for k in 1:ncol
        
        # residauls
        res.i <- dat.i - exp.i
        
        # permutation
        dist.i <- dist.obs
        p.i <- 1
        
        indiv.i <- 1:nrow(res.i)
        group1 <- c()
        group2 <- c()
        
        for (j in 1:NPerm) {
            # randoms sample and create data
            rand.j <- sample(indiv.i)
            dat.j <- dat.i[rand.j, ]
            
            # estimating parameters
            m.j <- lm(dat.j ~ group.i)
            
            # least-square mean
            obs.j <- by(m.j$fitted.value, group.i, colMeans)
            obs.j <- rbind(obs.j$"1", obs.j$"-1")
            
            # difference between groups
            obs.diff.j <- obs.j[1, ] - obs.j[2, ]
            dist.obs.j <- sqrt(t(obs.diff.j) %*% obs.diff.j)
            
            # record results
            dist.i <- c(dist.i, dist.obs.j)
            
            p.j <- ifelse(dist.obs.j >= dist.obs, 1, 0)
            
            p.i <- c(p.i, p.j)
            
            group1 <- c(group1, sum(obs.j[1, ]))
            group2 <- c(group2, sum(obs.j[2, ]))
            
        }  #end for j in NPerm
        results.i <- c(results.i, sum(p.i)/length(p.i))
        dists.i <- c(dists.i, dist.obs)
        
    }  #end for i in 1:nrow(comp)
    
    resultsP <- cbind(comp, results.i, dists.i)
    
    resultsP <- as.data.frame(resultsP)
    
    colnames(resultsP) <- c("Group 1", "Group 2", "Pr", "Distance")
    
    return(resultsP)
    
}  #end function permute.lda

