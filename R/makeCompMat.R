#' A function to create a pairwise comparison matrix
#' 
#' This function creates a pairwise comparison matrix for n groups.  All
#' possible pairwise combinations are created, with rows in the matrix equal to
#' the desired comparison.
#' 
#' 
#' @param ng A single number indicating the total number of unique groups
#' @return Returns a matrix with two columns and ng choose 2 rows.
#' @seealso \code{\link{PermuteLDA}}
#' @examples
#' 
#' makeCompMat(3)
#' 
#' makeCompMat(4)
#' 
#' data(Groups)
#' NGroups<-length(unique(Groups))
#' 
#' makeCompMat(NGroups)
#' 
#' @export makeCompMat
makeCompMat <- function(ng) {
    comparisons <- c()
    for (i in 1:ng) {
        c.i <- numeric(0)
        for (j in i:ng) {
            c.j <- c(i, j)
            if ((i - j) == 0) {
                next
            } else {
                c.i <- c(c.i, c.j)
            }  #end if/else i-j == 0
        }  #end for j in i:ng
        comparisons <- c(comparisons, c.i)
    }  #end for i in 1:ng
    
    compare <- matrix(comparisons, ncol = 2, byrow = TRUE)
    
    return(compare)
}  #end function makeCompMat

