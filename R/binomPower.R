binomPower<-function(ndads,mm,vv,tau2,nperms,nsims,nbins,doPlot=FALSE){
	getP<-function(ndads,mm,vv,tau2,nperms, nsims, nbins){
		ps<-c()
		for(i in 1:nsims){
			p.i<-simPower(ndads, mm, vv, tau2, nperms, nbins)
			ps<-c(ps,p.i)
		}
			return(ps)
		}#end getP
	
	#power analysis code
	params<-expand.grid(tau2,ndads,mm,vv)
	colnames(params)<-c('tau2','ndads','mm','vv')
	rm<-which(params$vv<params$mm)
	if(length(rm)>0){
		params<-params[-rm,]
	}
	tstamp<-as.character(as.integer(Sys.time()))
	write.csv(params,file=paste0('simPowerParams','_',tstamp,'.csv'),row.names=FALSE)
	
	ptest<-vector("list",nrow(params))
	out<-c()
	for(i in 1:nrow(params)){
		ptest.i<-getP(ndads=params$ndads[i], mm=params$mm[i], vv=params$vv[i], tau2=params$tau2[i], nperms=nperms,nsims=nsims,nbins=3)
		ptest[[i]]<-ptest.i
		out<-c(out,ptest.i)
	}
	
	out.save<-matrix(out,ncol=nsims,byrow=TRUE)
	write.csv(out.save,file=paste0(tstamp,'.csv'),row.names=FALSE)
	
	rejectM0<-function(x){
		return(length(which(x<0.05)))
	}
	
	reject<-apply(out.save,1,rejectM0)
	reject<-reject/nsims
	
	trueM0<-which(params$tau2[1:length(reject)]==0)
	trueM1<-(1:length(reject))[-trueM0]
	
	falseP<-rep(NA, length(reject))
	falseP[trueM0]<-reject[trueM0]
	falseP[trueM1]<-1-reject[trueM1]
	trueP<-1-falseP
	
	dat.plot<-data.frame(falseP,trueP, params$ndads[1:length(reject)], params$tau2[1:length(reject)])
	colnames(dat.plot)<-c('falseP','trueP','ndads','trueTau2')
	dat.plot$group<-paste0('trueTau2',':',dat.plot$trueTau,'-','nDads',':',dat.plot$ndads)
	dat.plot$trueTau<-as.factor(dat.plot$trueTau)
	dat.plot$ndads<-as.factor(dat.plot$ndads)
	write.csv(dat.plot,paste0('roc','_',tstamp,'.csv') ,row.names=FALSE)
	
	if(doPlot==TRUE){
		plotBinomPower(dat.plot,params,tstamp)
	}
	
	out<-list('roc'=dat.plot,'params'=params,'results'=out.save)
	
	return(out)
}