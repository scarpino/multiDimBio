makeCompMat <-
function(ng){
	comparisons<-c()
	for(i in 1:ng){
		c.i<-numeric(0)
		for(j in i:ng){	
		c.j<-c(i,j)
		if((i-j)==0){
			next}else{
		c.i<-c(c.i,c.j)
				}#end if/else i-j == 0
			}#end for j in i:ng
		comparisons<-c(comparisons,c.i)
		}#end for i in 1:ng
		
	compare<-matrix(comparisons,ncol=2,byrow=TRUE)
	
	return(compare)
	}#end function makeCompMat

