Power <-
function(func='PermuteLDA', N='DEFAULT.N', effect.size='DEFAULT.e',trials=100){
	cat('Power Analysis -',func,'\n','\n')
	if(length(N) <2) {stop('must have more than two group sizes')}
	if(func =='FSelect' & min(N <6)){stop('must have more than 5 indiv with FSelect')}
	if(is.character(effect.size)==TRUE){
		ES<-c(0.1,0.4,0.8,1.6)
		}else{
			ES<-effect.size
			}
	if(is.character(N)==TRUE){
		N<-c(6,12,24,48,96)
		}
	
	counter<-length(ES)
		
	for(e in ES){
		start<-Sys.time()
		
		effect<-e
		
		RESULTS<-c()
		
		for(n in N){
			RESULTS.P<-c()
			
			for(k in 0:4){
				r.p<-c()
					for(i in 1:trials){
						if(k==0){
							mu<-c(0,0,0,0)}else{
								if(k==1){
									mu<-c(effect,0,0,0)}else{
									if(k==2){
									mu<-c(effect,effect,0,0)}else{
										if(k==3){
									mu<-c(effect,effect,effect,0)}else{mu<-c(effect,effect,effect,effect)
										}#end if/else (k=3)
									}#end if/else (k=2)
								}#end if/else (k=1)
							}#end if/else (k=0)
						t1<-rnorm(n)
						t2<-c(rnorm(n/2,mu[4]),rnorm(n/2,0))
						t3<-c(rnorm(n/2,mu[1]),rnorm(n/2,0))
						t4<-c(rnorm(n/2,mu[2]),rnorm(n/2,0))
						t5<-c(rnorm(n/2,mu[3]),rnorm(n/2,0))
						t6<-runif(n)

						#data
						T<-cbind(t1,t2,t3,t4,t5,t6)

						#groups
						GR<-rep(1:2,each=n/2)

						#results
						if(func=='PermuteLDA'){
							func.i<-PermuteLDA(T,GR,100)
							if(func.i[3]<=0.05){
								p.i<-1}else{
								p.i<-0
								}#end if/else(func.i[3])
						}else{
							func.i<-FSelect(T,GR,4)
							ks<-	0:k
							if(length(k)==0){
								is.sig<-sum(func.i$PrF<=0.05)
								}else{
								ks<-ks[-1]
								is.sig<-sum(func.i$PrF[ks]<=0.05)
							}#end if/else length(k)
							
							if(is.sig>=1){
								p.i<-1}else{
									p.i<-0
								}#end if/else is.sig>=1
						}#end if/else func==
							
							r.p<-c(r.p,p.i)
						}#end for i

				RESULTS.P<-c(RESULTS.P,sum(r.p)/i)
			}#end for k

		RESULTS<-c(RESULTS,RESULTS.P)
	}#end for n

	RP<-matrix(RESULTS,ncol=5,byrow=TRUE)

	X<-N
	X<-X/2
	Size<-rep(X, times=5)
	SigTraits<-c('0','1','2','3','4')
	SigTraits<-rep(SigTraits,each=length(X))
	RES.P<-matrix(RP, ncol=1)
	DAT<-data.frame(Size,RES.P,SigTraits)
	
	pal<-palette(brewer.pal(5,"YlGnBu"))
	pal<-palette(brewer.pal(5,"YlGnBu"))
	pal[1]<-'black'
	
	p<-ggplot(DAT,aes(Size,RES.P))
	p.save<-p+geom_line(aes(group=SigTraits,colour=SigTraits),size=2)+scale_x_continuous(expand=c(0,0),name='Individuals Per Group')+scale_y_continuous(expand=c(0,0),limits=c(0,1),name='Proportion Significant')+geom_hline(aes(yintercept=0.8),linetype=3)+scale_colour_manual(values=pal)+theme(legend.position ='right', legend.background=element_rect(fill='#ffffffaa',colour='black'),panel.background=element_rect(fill='white',colour='black'),axis.text=element_text(colour='black',size=15),axis.title=element_text(colour='black',size=15), legend.key=element_rect(fill='white '),panel.grid.minor=element_blank(),panel.grid.major=element_blank())

	ggsave(filename=paste(as.character(effect),func,'pdf',sep='.'), p.save)
	
	write.csv(DAT,file=paste(as.character(effect),func,'csv',sep='.'))
	
	counter<-counter-1
	end<-Sys.time()
	total<-end-start
	remaining<-total*counter
	print(remaining)
	cat('The above was the estimated remaining time','\n','\n')
}#end for e
cat('Results saved as .pdfs in working dir','\n','\n')
}#end FUNCTION

