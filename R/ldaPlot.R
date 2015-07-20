ldaPlot <-
function(Data, Groups, palette='BrBG', axes=c(1,2,2,3,1,3)){
	LD1<-NULL
	LD2<-NULL
	
	if(min(Data, na.rm=TRUE)<1e-4){
	TOL<-min(Data, na.rm=TRUE)
	}else{
		TOL<-1e-4
		}
			
	Groups<-as.factor(Groups)
	Groups<-as.numeric(Groups)
	
	if(length(unique(Groups))>2){
		
	if(length(unique(Groups))>6){
	cat('Will not work with more than 6 groups','\n','\n')
	}#end if >6

	ld1<-lda(Groups~Data, tol=TOL)
	pm1<-predict(ld1)
	GC<-as.factor(pm1$class)
	pm1df<-data.frame(GC,pm1$x)
	
	cat(paste('Plotting specified axes using',palette,'palette',sep=' '),'\n','\n')
	
	AX<-matrix(axes,ncol=2,byrow=TRUE)
	
	for(i in 1:nrow(AX)){
	AX1<-AX[i,1]+1
	AX2<-AX[i,2]+1
	pm1df.i<-pm1df[,c(1,AX1,AX2)]
	colnames(pm1df.i)<-c('GC','LD1','LD2')
	#pi<-ggplot(pm1df,aes(pm1df[,AX1],pm1df[,AX2]))
	pi<-ggplot(pm1df.i,aes(LD1,LD2))
	pi.save<-pi+geom_point(aes(colour=GC, shape=GC), size=3.5)+scale_colour_brewer(palette=palette)+scale_x_continuous(name=colnames(pm1df)[AX1])+scale_y_continuous(name=colnames(pm1df)[AX2])+theme(legend.position ='right', legend.background=element_rect(fill='#ffffffaa',colour='black'),panel.background=element_rect(fill='white',colour='black'),axis.text=element_text(colour='black',size=15),axis.title=element_text(colour='black',size=15), legend.key=element_rect(fill='white '),panel.grid.minor=element_blank(),panel.grid.major=element_blank())
	
	timestamp<-as.character(as.integer(Sys.time()))
	ggsave(paste(timestamp,'ldaPlot-AllGroups.pdf'),pi.save)
	}#end for i
	}#end if unique(groups)>2
	
	#1 comparison matrix
	ngroups<-length(unique(Groups))
	comp<-makeCompMat(ngroups)
	
	for(j in 1:nrow(comp)){
		use.j<-which(Groups==comp[j,1]|Groups==comp[j,2])
		
		ld.j<-lda(Groups[use.j]~Data[use.j,], tol=TOL)
		pm.j<-predict(ld.j)
		GC.j<-as.factor(pm.j$class)
		pmdf.j<-data.frame(GC.j,pm.j$x)
		colnames(pmdf.j)<-c('GC','LD1')	
		
		p.j<-ggplot(pmdf.j, aes(x=LD1))
		
		psave.j<-p.j+geom_density(aes(fill=GC),xlab='Linear Discriminant 1 Score')+scale_fill_manual(values=c('#8C510A66','#01665E66'))+geom_point(aes(y=0,x=LD1,shape=GC),size=3.5)+geom_point(aes(y=0,x=LD1,colour=GC,shape=GC))+scale_colour_manual(values=c('#8C510A','#01665E'))+xlab('Linear Discriminant 1 Score')+ylab('Density')+theme(legend.position ='right', legend.background=element_rect(fill='#ffffffaa',colour='black'),panel.background=element_rect(fill='white',colour='black'),axis.text=element_text(colour='black',size=15),axis.title=element_text(colour='black',size=15), legend.key=element_rect(fill='white '),panel.grid.minor=element_blank(),panel.grid.major=element_blank())
		
		timestamp<-as.character(as.integer(Sys.time()))
		ggsave(paste(timestamp,comp[j,1],comp[j,2],'Group LDA Plot.pdf',sep=' '), psave.j)
		}#end for j
		
}#end FUNCTION

