#' Estimates the heritability of a binomial trait
#' 
#' Estimates the narrow-sense heritability of a binomial trait and calculates a
#' p value by randomization.
#' 
#' Estimates the narrow-sense heritability of a binomial trait.  This function
#' works by fitting two models, one with and one without a random-effect of
#' sire.  These models are compared by randomizing the sire ids nreps times and
#' re-fitting the model.  For each of the nreps model pairs, a deviance is
#' calculated and a "p value" estimated by comparing that distribution of
#' deviance to the observed.  The heritability is approximatly
#' tau2/(tau2+(pi/sqrt(3))^2), where tau2 is the random-effect variance due to
#' sire.
#' 
#' @param data a (non-empty) numeric matrix with three columns.  The first two
#' should contain the trait data (number of occurances of each outcome type)
#' and the third should contain the group ids.
#' @param nreps a (non-empty) numeric value indicating the number of resamples
#' to perform when calculating the emperical p value.
#' @return Returns a list.  The list contains: \item{h2}{ The estimated
#' narrow-sense heritability.  The narrow-sense heritability is approximatly
#' tau2/(tau2+(pi/sqrt(3))^2), where tau2 is the random-effect variance due to
#' sire. }
#' 
#' \item{pval}{ The probability that the best-fit model includes an extra
#' variance term for sire (random effect of dad).  The value is calculated by
#' comparing the deviances from nreps number of randomized model comparisions.
#' }
#' 
#' \item{deviance}{ The deviance between a null model without a random effect
#' of dad and a model with. }
#' 
#' \item{sim}{ The simulated deviances used in calculating the p value in pval.
#' }
#' 
#' \item{obsMod}{ The glmer model object resulting from the observed data. }
#' @examples
#' 
#' 	ndads <- 18
#' 	mm <- 4
#' 	vv <-  6
#' 	tau2 <- 1.5
#' 	nbins <- 3
#' 	
#' 	mylogit <- function(x) log(x/{1-x})
#' 	ilogit <- function(x) 1/{1+exp(-x)}
#' 	
#' 	swimprob <- ilogit(rnorm(ndads, 0, sqrt(tau2)))
#' 	mytable <- NULL
#' 	for(i in 1:ndads) {
#' 		bincounts <- pmax(1,rnbinom(nbins, mu = mm, size = mm^2/{vv-mm}))
#' 		swim <- rbinom(3, bincounts,swimprob[i])
#' 		set <- bincounts - swim
#' 		theserows <- data.frame(set=set,swim=swim, Dad = i, Bin = 1:nbins)
#' 		mytable <- rbind(mytable, theserows)
#' 	}
#' 	
#' 	est <- h2Estimate(mytable,nreps=10)
#' 	
#' 	print(est$h2)
#' 
#' @export h2Estimate
h2Estimate <- function(data, nreps = 1000) {
    colnames(data) <- c("trait1", "trait2", "dad")
    hm0.real = glm(cbind(trait1, trait2) ~ 1, data = data, family = binomial)
    # hm1.real = glmer(cbind(trait1, trait2) ~ (1 | dad),
    # data=data, family=binomial, REML=FALSE)
    hm1.real = glmer(cbind(trait1, trait2) ~ (1 | dad), data = data, 
        family = binomial)
    Dobs = as.numeric(deviance(hm0.real) - deviance(hm1.real))
    
    ncases <- nrow(data)
    perm1 = replicate(nreps, {
        neworder = sample(1:ncases, ncases)
        ptable = data.frame(dad = data$dad, t1 = data$trait1[neworder], 
            t2 = data$trait2[neworder])
        hm0 = glm(cbind(t1, t2) ~ 1, data = ptable, family = binomial)
        # hm1 = glmer(cbind(t1,t2) ~ (1 | dad), data=ptable,
        # family=binomial, REML=FALSE)
        hm1 = glmer(cbind(t1, t2) ~ (1 | dad), data = ptable, 
            family = binomial)
        D = as.numeric(deviance(hm0) - deviance(hm1))
        D
    })
    
    pval <- length(which(perm1 > Dobs))/length(perm1)
    
    trueTau2 <- (VarCorr(hm1.real)$dad[1])^2
    h2 <- 4 * (trueTau2/(trueTau2 + (pi/sqrt(3))^2))
    
    out <- list(h2 = h2, pval = pval, deviance = Dobs, sim = perm1, 
        obsMod = hm1.real)
    
    return(out)
}
