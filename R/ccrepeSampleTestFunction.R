ccrepeSampleTestFunction  <- function
(x, y = NA){
#*************************************************************************************
#*  	ccrepeSampleTestFunction                                                     *
#* This a simple example fot a test measurent function to be used with ccrepe        *
#* used in the same fashion that cor would be used                                   *
#* Some properties of the function:                                                  *
#*                                                                                   *
#* 1. Be able to take either two inputs which are vectors or one input which         *
#*   is either a matrix or a data frame                                              *
#* 2.In the case of two inputs, return a single number                               *
#* 3.In the case of one input, return a matrix in which the (i,j)th entry            *
#*   is the similarity score for column i and column j in the original matrix        *
#* 4.Resulting matrix must be symmetric                                              *
#* 5.The inputs must be named x and y                                                *
#*************************************************************************************
	if(is.vector(x) && is.vector(y)) return(.5)                                   
	if(is.matrix(x) && is.na(y)) return(matrix(rep(.5,ncol(x)^2),ncol=ncol(x)))   
	if(is.data.frame(x) && is.na(y)) return(matrix(rep(.5,ncol(x)^2),ncol=ncol(x)))
	else stop('ERROR')
	}