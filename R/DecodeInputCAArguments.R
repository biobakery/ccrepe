DecodeInputCAArguments <-
function(CA){
#*************************************************************************************
#*	Decode and validate input parameters                                             *
#*************************************************************************************
	if (length(CA$data1) == 1)					#If the user did not select at least one input dataframe - Stop the run
		{
		if (is.na(CA$data1))
			stop('You must select at least one input data frame (data1<-YourData)') 
		}
	CA$OneDataset <- FALSE						#Set the symmmetrical flag to False (Will be true if data2=data1)		
	if (length(CA$data2) == 1)					#If user did not enter data2 - we assume he wants data2=data1
		if (is.na(CA$data2))					
		{
		CA$OneDataset <- TRUE					#And set up the Symmetrical flag
		}


	if ( CA$subset.cols.1[1]  == 0)					#If user did not pass any request for the number fo cols in data1 to process
		CA$subset.cols.1<-c(1:ncol(CA$data1))			# - We will use All the columns	

	

	if (CA$OneDataset == FALSE)
		{
		if ( CA$subset.cols.2[1]  == 0)					#If user did not pass any request for the number fo cols in data1 to process
		CA$subset.cols.2<-c(1:ncol(CA$data2))				# - We will use All the columns	
		}
	
	
	MethodsTable <- c('spearman','kendall','pearson')
	if   (CA$method  %in% MethodsTable == FALSE)
		{
		cat('\nInvalid correlation method selected : ',CA$method,' - using pearson instead\n')
		CA$method <- 'pearson'
		}

	if    (!is.na(CA$outdist))							#If the user passed a file - open it
		{
		CA$outdistFile = file(CA$outdist,open='at')		#Open outdist file
		}
	return(CA)			 				#Return list of decoded input parameters
}
