simPower <- function(ndads, mm, vv, tau2, nperms, nbins) {
    mylogit = function(x) log(x/{
        1 - x
    })
    ilogit = function(x) 1/{
        1 + exp(-x)
    }
    
    swimprob = ilogit(rnorm(ndads, 0, sqrt(tau2)))
    mytable = NULL
    for (i in 1:ndads) {
        bincounts = pmax(1, rnbinom(nbins, mu = mm, size = mm^2/{
            vv - mm
        }))
        swim = rbinom(3, bincounts, swimprob[i])
        set = bincounts - swim
        theserows = data.frame(Dad = i, Bin = 1:nbins, set = set, 
            swim = swim)
        mytable = rbind(mytable, theserows)
    }
    ncases = nrow(mytable)
    empfreq = aggregate(swim ~ Dad, data = mytable, sum)[, 2]/aggregate(swim + 
        set ~ Dad, data = mytable, sum)[, 2]
    
    # Compute the likelihood ratio statistic
    hm0 = glm(cbind(swim, set) ~ 1, data = mytable, family = binomial)
    hm1 = glmer(cbind(swim, set) ~ (1 | Dad), data = mytable, 
        family = binomial, REML = FALSE)
    Dsim = as.numeric(deviance(hm0) - deviance(hm1))
    
    perm1 <- c()
    for (i in 1:nperms) {
        neworder <- sample(1:ncases, ncases)
        ptable <- data.frame(Dad = mytable$Dad, Bin = mytable$Bin, 
            swim = mytable$swim[neworder], set = mytable$set[neworder])
        hm0 <- glm(cbind(swim, set) ~ 1, data = ptable, family = binomial)
        hm1 <- glmer(cbind(swim, set) ~ (1 | Dad), data = ptable, 
            family = binomial)
        D <- as.numeric(deviance(hm0) - deviance(hm1))
        perm1 <- c(perm1, D)
    }
    
    pval <- length(which(perm1 > Dsim))/length(perm1)
    return(pval)
}
