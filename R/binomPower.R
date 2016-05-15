#' Power analysis for estimating the heritability of a binomial trait
#' 
#' Performs a power analysis for estimating the heritability of a binomial
#' trait.  This function can take a long time to run if either nsims or nperms
#' is large.
#' 
#' 
#' @param ndads a (non-empty) numeric value indicating the number of dads.
#' @param mm a (non-empty) numeric value indicating the mean number of
#' offspring per dad per bin (normal dist).
#' @param vv a (non-empty) numeric value indicating the variance in offspring
#' per dad per bin (normal dist).
#' @param tau2 a (non-empty) numeric value indicating the dad effect
#' (narrow-sense heritability ~ tau2/(tau2+(pi/sqrt(3))^2)).
#' @param nperms a (non-empty) numeric value indicating the number of bootstrap
#' permutations to use for caluclating a p value.
#' @param nsims a (non-empty) numeric value indicating the number of
#' simulations to run per parameter combination.
#' @param nbins a (non-empty) numeric value indicating the number of bins, data
#' are pooled before analysis.
#' @param doPlot a (non-empty) logical value indicating whether to plot the
#' results of the power analysis.
#' @return Returns a list and an optional set of .pdfs (if doPlot==TRUE).  The
#' list contains: \item{roc}{ a data.frame with the summarized results of the
#' power analysis. }
#' 
#' \item{params}{ a numeric matrix with the paramater values. }
#' 
#' \item{results}{ a numeric matrix with the full results of the analysis. }
#' @examples
#' 
#' 	ndads <- c(9,18)
#' 	mm <- 4.629634
#' 	vv <- 6.31339
#' 	tau2 <- c(0,0.5)
#' 	nperms <- 2
#' 	nsims <- 2
#' 	nbins <- 3
#' 	doPlot <- TRUE
#' 	binomPower(ndads,mm,vv,tau2,nperms,nsims,nbins,doPlot)
#' 
#' @export binomPower
binomPower <- function(ndads, mm, vv, tau2, nperms, nsims, nbins, alpha = 0.05, 
    doPlot = FALSE) {
    
    # power analysis code
    params <- expand.grid(tau2, ndads, mm, vv)
    colnames(params) <- c("tau2", "ndads", "mm", "vv")
    rm <- which(params$vv < params$mm)
    if (length(rm) > 0) {
        params <- params[-rm, ]
    }
    
    ptest <- vector("list", nrow(params))
    out <- c()
    for (i in 1:nrow(params)) {
        ptest.i <- getP(ndads = params$ndads[i], mm = params$mm[i], 
            vv = params$vv[i], tau2 = params$tau2[i], nperms = nperms, 
            nsims = nsims, nbins = 3)
        ptest[[i]] <- ptest.i
        out <- c(out, ptest.i)
    }
    
    out.save <- matrix(out, ncol = nsims, byrow = TRUE)
    
    reject <- apply(out.save, 1, function(x, alpha) return(length(which(x < alpha))), alpha = alpha)
    reject <- reject/nsims
    
    trueM0 <- which(params$tau2[1:length(reject)] == 0)
    trueM1 <- (1:length(reject))[-trueM0]
    
    falseP <- rep(NA, length(reject))
    falseP[trueM0] <- reject[trueM0]
    falseP[trueM1] <- 1 - reject[trueM1]
    trueP <- 1 - falseP
    
    dat.plot <- data.frame(falseP, trueP, params$ndads[1:length(reject)], 
        params$tau2[1:length(reject)])
    colnames(dat.plot) <- c("falseP", "trueP", "ndads", "trueTau2")
    dat.plot$group <- paste0("trueTau2", ":", dat.plot$trueTau, 
        "-", "nDads", ":", dat.plot$ndads)
    dat.plot$trueTau <- as.factor(dat.plot$trueTau)
    dat.plot$ndads <- as.factor(dat.plot$ndads)
    
    tstamp <- as.numeric(Sys.time())
    if (doPlot == TRUE) {
        plotBinomPower(dat.plot, params, tstamp)
    }
    
    out <- list(roc = dat.plot, params = params, results = out.save)
    
    return(out)
}
