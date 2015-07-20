ZTrans <-
function(DATA){
	cmean<-apply(DATA,2,mean,na.rm=TRUE)
	cmat<-matrix(cmean,ncol=ncol(DATA),nrow=nrow(DATA),byrow=TRUE)
	muDAT<-DATA-cmat
	
	csd<-apply(DATA,2,sd,na.rm=TRUE)
	cdmat<-matrix(csd,ncol=ncol(DATA),nrow=nrow(DATA),byrow=TRUE)
	
	Z.SCORES<-muDAT/cdmat
	
	return(Z.SCORES)
	}#end FUNCTION

