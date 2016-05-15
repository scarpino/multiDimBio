PermuteLDA <- function(Data, Groups, NPerm, Missing.Data = "Complete") {
    
    Groups <- as.factor(Groups)
    Groups <- as.numeric(Groups)
    
    # Missing data
    if (Missing.Data == "Complete" & length(which(is.na(Data) == 
        TRUE)) > 0) {
        cat("Missing data imputed using CompleteData to exlude missing data set Missing.Data=Remove", 
            "\n", "\n")
        Data <- CompleteData(Data, NPCS = floor(ncol(Data)/5), 
            cut.trait = 1, cut.ind = 1)
    }
    
    if (Missing.Data == "Remove" & length(which(is.na(Data) == 
        TRUE)) > 0) {
        cat("Individuals with missing data removed, to impute missing data set Missing.Data=Complete", 
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

