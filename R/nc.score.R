#!/usr/bin/env Rscript
nc.score <-
function(
#*************************************************************************************
#*  	nc.score                                                                     *
#*************************************************************************************
   	x=NA,						#First input 
	y=NA,						#Second input
	adBins=NA,						#Number of Input Bins
	min.abundance=NA,				#Minimum Abundance
	min.samples=NA)				#Minimum Samples
	
{
#*************************************************************************************
#*  	NC_Score                                                                     *
#*************************************************************************************
	if (is.numeric(x)		    #Is it a vector? 		
		& length(x) > 1
		& is.numeric(y) 
		& length(y) > 1 
		& length(x) == length(y)
		& class(x) == "numeric"
	    & class(y) == "numeric"
		)
		{
			x.discretized = as.matrix(discretize(x,nbins = sqrt(length(x))))	#Discretize 
			y.discretized = as.matrix(discretize(y,nbins = sqrt(length(y))))	#Discretize
			nc.score.result = nc.score.helper(x.discretized,y.discretized)					#Invoke the function
			nc.score.result = nc.score.renormalize (x.discretized, y.discretized, nc.score.result)  #Normalize the results 
			return(nc.score.result)
		}
		
	#************************************************************************
	#*        It is a dataframe or a matrix                                 *
	#************************************************************************

	CA = preprocess_nc_score_input (
			x, 										#First Input
			y,										#Second Input
			adBins,									#Bins
			min.abundance,							#Minimum Abundance
			min.samples)							#Minimum Samples
	
	x <- CA$x										#Get it from common area
	y <- CA$y										#Get it from common area
    x.discretized <- CA$x.discretized				#Get it from Common Area
	y.discretized <- CA$y.discretized				#Get it from Common Area

	for (i in 1:ncol(x))												#Loop on the columns of the first matrix
		{
			x.discretized[,i] = discretize(x[,i],nbins=sqrt(length(x[,i])))[,1]	#Discretize it and post the value in the discretized matrix
		}
	for (i in 1:ncol(y))												#Do the same for the second matrix
		{
			y.discretized[,i] = discretize(y[,i],nbins=sqrt(length(y[,i])))[,1]	#Post it in the second discretized matrix
		}
	nc.score.result <- matrix(nrow=ncol(x.discretized),ncol=ncol(y.discretized))		#Build results area
	
	for (i in 1:ncol(x.discretized))									#Loop on the columns of the first matrix
		{
			for (j in 1:ncol(y.discretized))							#Loop on the columns of the second matrix
				{
					nc.score.result[i,j]<-nc.score.helper(x.discretized[,i],y.discretized[,j])	#Post it in the results area
				}			
		}
 	nc.score.result <- nc.score.renormalize (x.discretized, y.discretized, nc.score.result)  #Normalize the results 
								
	
return(nc.score.result)
}