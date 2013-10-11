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
			nonmissing.subjects <- intersect(which(!is.na(x)),which(!is.na(y)))      #Find subjects which are present in both x and y
			CA$x <- x[nonmissing.subjects]					#Post x without missing subjects to common area
			CA$y <- y[nonmissing.subjects]					#Post y without missing subjects to common area
			CA$verbose <- verbose						#Post the verbose flag
 
			CA <- process.input.bins(bins, CA)			#Process bins input request
			if (length(CA$bins) == 1)					#Check if bins is a number or a vector with entries
				{
					x.discretized = as.matrix(discretize(CA$x,nbins = CA$bins))	#Discretize x
					y.discretized = as.matrix(discretize(CA$y,nbins = CA$bins))	#Discretize y
				} else
				{
					x.discretized  = as.matrix(findInterval(CA$x,CA$bins))	 #Use the bins provided by the User
					y.discretized  = as.matrix(findInterval(CA$y,CA$bins))	 #Use the bins provided by the User
				}
			

			nc.score.result = nc.score.vectors.helper(x.discretized,y.discretized,CA)		#Invoke the function
			nc.score.result = nc.score.renormalize (x.discretized, y.discretized, nc.score.result)  #Normalize the results 
			CA$nc.score.result <- nc.score.result		#Post the result
			if (CA$verbose == TRUE)
				{return(CA)}
				else
				{return(nc.score.result)}	
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
					x.discretized[,i] = findInterval(x[,i],CA$bins)	 #Use the bins provided by the User
				}
		}

 	CA$nc.score.matrix <-nc.score.helper(x.discretized,CA)	#Post it in the results area
 	CA$nc.score.matrix <- nc.score.renormalize (x.discretized, NA, CA$nc.score.matrix)  #Normalize the results 
	rownames(CA$nc.score.matrix)  <- colnames(CA$x)	#Post the row names
	colnames(CA$nc.score.matrix)  <- colnames(CA$x)	#Post the col names
	CA$x.discretized <- NULL		#Not needed anymore
	CA$x <- NULL					#Not Needed anymore
	
	diag(CA$nc.score.matrix)<-NA	#We are setting the diagonal entries in the matrix to NA
	if (length(CA$columns.not.passing.qc) > 0)  #If there were columns that did not pass QA, we need to add corr cols with NA
		{
			original.nc.score.dim <- ncol(CA$nc.score.matrix)  #Columns in the original matrix
			na.c <-  rep(NA, original.nc.score.dim)   	#A vector of NAs to be inserted
			rebuilt.matrix <- CA$nc.score.matrix					#Allocate the rebuilt matrix
			for (indx in 1:length(CA$columns.not.passing.qc))
				{
					left.part <- rebuilt.matrix[,1:CA$columns.not.passing.qc [indx] -1]		#Left part of the rebuilt matrix
					right.part <- rebuilt.matrix[,CA$columns.not.passing.qc [indx]:ncol(rebuilt.matrix)]		#Right part of the rebuilt matrix
					rebuilt.matrix <- cbind(left.part,na.c,right.part)		#Insert column of NAs
				}

			na.r <- 	rep(NA, ncol(rebuilt.matrix))   	#A vector of NAs to be inserted as a row
			for (indx in 1:length(CA$columns.not.passing.qc))
				{
				upper.part <- rebuilt.matrix[1:CA$columns.not.passing.qc [indx] -1,]		#Upper part of the rebuilt matrix
				lower.part <- rebuilt.matrix[CA$columns.not.passing.qc [indx]:nrow(rebuilt.matrix),]		#Right part of the rebuilt matrix
				rebuilt.matrix <- rbind(upper.part,na.r,lower.part)		#Insert row of NAs
				}
			CA$nc.score.matrix <- rebuilt.matrix				#Post the matrix
		}
 
	CA$input.total.cols <- NULL		#Not needed anymore
	CA$columns.not.passing.qc <- NULL  #Not needed anymore
	
	if (!CA$verbose == TRUE)		#If abbreviated output
		{
		CA <- CA$nc.score.matrix	#just post the resulting matrix
		}

return(CA)
}