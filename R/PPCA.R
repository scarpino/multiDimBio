PPCA<-
function(Data,nPCs=4,CENTER=TRUE,SCALE='vector'){
	PC<-pca(Data,nPcs=nPCs, method='ppca', center=CENTER,scale=SCALE)
	return(PC)
	}#end Function
	