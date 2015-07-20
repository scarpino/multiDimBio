h2Estimate<-function(data,nreps=1000){
	colnames(data)<-c('trait1','trait2','dad')
	hm0.real = glm(cbind(trait1, trait2) ~ 1, data=data, family=binomial)
	#hm1.real = glmer(cbind(trait1, trait2) ~ (1 | dad), data=data, family=binomial, REML=FALSE)
	hm1.real = glmer(cbind(trait1, trait2) ~ (1 | dad), data=data, family=binomial)
	Dobs = as.numeric(deviance(hm0.real) - deviance(hm1.real))
	
	ncases<-nrow(data)
	perm1 = replicate(nreps, {
		neworder = sample(1:ncases, ncases)
		ptable = data.frame(dad = data$dad, t1 = data$trait1[neworder], t2 = data$trait2[neworder])
		hm0 = glm(cbind(t1,t2) ~ 1, data=ptable, family=binomial)
		#hm1 = glmer(cbind(t1,t2) ~ (1 | dad), data=ptable, family=binomial, REML=FALSE)
		hm1 = glmer(cbind(t1,t2) ~ (1 | dad), data=ptable, family=binomial)
		D = as.numeric(deviance(hm0) - deviance(hm1));
	D;
	}
	)

	pval<-length(which(perm1>Dobs))/length(perm1)
	
	trueTau2<-(VarCorr(hm1.real)$dad[1])^2
	h2<- 4 * (trueTau2/(trueTau2 +(pi/sqrt(3))^2))
	
	out<-list('h2'=h2, 'pval'=pval, 'deviance'=Dobs, 'sim'=perm1, 'obsMod'=hm1.real)
	
	return(out)
}