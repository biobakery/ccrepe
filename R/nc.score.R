#!/usr/bin/env Rscript
nc.score <-
function(
#*************************************************************************************
#*  	nc.score                                                                     *
#*************************************************************************************
   	x=NA,						#First input 
	y=NA,						#Second input
	bins=NA,					#Number of Input Bins
	verbose = FALSE,			#Request for verbose output?
	min.abundance=0.0001,		#Minimum Abundance
	min.samples=0.1)			#Minimum Samples
	
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
			CA <-list()									#Set the common area
			CA$x <- x									#Post x to common area
			CA$verbose <- verbose						#Post the verbose flag
			CA <- process.input.bins(bins, CA)
			x.discretized = as.matrix(discretize(x,nbins = CA$bins))	#Discretize x
			y.discretized = as.matrix(discretize(y,nbins = CA$bins))	#Discretize y
			nc.score.result = nc.score.vectors.helper(x.discretized,y.discretized)					#Invoke the function
			nc.score.result = nc.score.renormalize (x.discretized, y.discretized, nc.score.result)  #Normalize the results 
			return(nc.score.result)
		}
		
	#************************************************************************
	#*        It is a dataframe or a matrix                                 *
	#************************************************************************
	CA = preprocess_nc_score_input (
			x, 										#First Input
			bins,									#Bins
			min.abundance,							#Minimum Abundance
			min.samples)							#Minimum Samples
	
	x <- CA$x										#Get the filtered x from common area
	CA$verbose <- verbose							#Post the verbose flag

    x.discretized <- CA$x.discretized				#Get it from Common Area
 
	for (i in 1:ncol(x))												#Loop on the columns of the first matrix
		{
			if (length(CA$bins) == 1)					#Check if bins is a number or a vector with entries
				{
					x.discretized[,i] = discretize(x[,i],nbins=CA$bins)[,1]	#Bins is a number; Discretize it and post the value in the discretized matrix
				}
				else
				{
					x.discretized[,i] = findInterval(x[,i],CA$bins)	+ 1	#Use the bins provided by the User
				}
		}

 	CA$nc.score.matrix <-nc.score.helper(x.discretized,CA)	#Post it in the results area
 	CA$nc.score.matrix <- nc.score.renormalize (x.discretized, NA, CA$nc.score.matrix)  #Normalize the results 
	rownames(CA$nc.score.matrix)  <- colnames(CA$x)	#Post the row names
	colnames(CA$nc.score.matrix)  <- colnames(CA$x)	#Post the col names
	CA$x.discretized <- NULL		#Not needed anymore
	CA$x <- NULL					#Not Needed anymore
	diag(CA$nc.score.matrix)<-NA	#We are setting the diagonal entries in the matrix to NA
	if (!CA$verbose == TRUE)		#If abbreviated output
		{
		CA <- CA$nc.score.matrix	#just post the resulting matrix
		}
return(CA)
}