BoxWhisker <-
function(DATA,GROUPS,palette='Paired'){
	
	message<-paste('Using',palette,sep=' ')
	
	cat(paste(message,'palette from colorbrewer2.org',sep=' '),'\n','\n')

	DAT<-as.matrix(DATA)
	SCORE<-matrix(DAT,ncol=1)
	TRAIT<-colnames(DATA)
	TRAIT<-rep(TRAIT,each=nrow(DAT))
	GROUP<-rep(GROUPS,ncol(DAT))
	
	DAT.df<-data.frame(SCORE,TRAIT,GROUP)
	DAT.df$GROUP<-factor(DAT.df$GROUP)
	DAT.TRAIT<-factor(DAT.df$TRAIT)
	
	p<-ggplot(DAT.df, aes(TRAIT,SCORE))
	p.save<-p+geom_boxplot(aes(fill=GROUP),outlier.shape=NA)+scale_fill_brewer(palette=palette)+theme(legend.position ='right', legend.background=element_rect(fill='#ffffffaa',colour='black'),panel.background=element_rect(fill='white',colour='black'),axis.text.y=element_text(colour='black',size=15),axis.text.x=element_text(colour='black',size=8),axis.title=element_text(colour='black',size=15), legend.key=element_rect(fill='white '),panel.grid.minor=element_blank(),panel.grid.major=element_blank())
	
	timestamp<-as.character(as.integer(Sys.time()))
	ggsave(filename=paste(timestamp,'BoxWhisker.pdf',sep='_'), p.save)
	}#end FUNCTION

