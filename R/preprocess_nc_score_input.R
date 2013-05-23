preprocess_nc_score_input <-
function(
#********************************************************************************************
#*  	Preprocess input and build the common area                                          *
#********************************************************************************************
			x, 										#First Input
			input.bins,
			input.min.abundance,
			input.min.samples)
{

	CA <-list()												#Common area definition 	 
  	if (!suppressWarnings(!is.na(as.numeric(input.min.abundance))))	#check if the User entered a valid threshold1
		{
		cat('\nInvalid min.abundance entered - using default=0.0001\n')
		CA$min.abundance = 0.0001
		 }				#If it is not valid - force default
	else
		{CA$min.abundance = input.min.abundance}	#Else - use it

	if (!suppressWarnings(!is.na(as.numeric(input.min.samples))))	#check if the User entered a valid threshold1
		{
		cat('\nInvalid min.samples entered - using default=0.1\n')
		CA$min.samples = 0.1 				#If it is not valid - force default
		}
	else
		{CA$min.samples = input.min.samples}	#Else - use it
 
 

	#*********************************************************
	#* Filter                                                *
	#*********************************************************
	x <-qc_filter(x,CA)										#filter by abundance throsholds
	x <-na.omit(x)											#remove NAs
	if (is.null(rownames(x))) {rownames(x)<-seq(1:nrow(x))} #If there are no row names - plug them in
	
	CA$x <-  x 												#Post it to common Area
	cat('\nAfter filtering there are ',nrow(CA$x),'samples with ',ncol(CA$x),
		'variables (bugs): Output is a matrix of ',ncol(CA$x), ' by ' , ncol(CA$x), '\n')
	CA$x.discretized <- matrix(nrow=nrow(CA$x),ncol=ncol(CA$x))			#Build x Discretized empty Matrix			
	#****************************************************************************************
	#*   Process Bins  passed by the user  or set default                                   *
	#****************************************************************************************
	CA <- process.input.bins(input.bins, CA)			#Process Input number of bins or set default
	
	return(CA)
}



process.input.bins <-
function(input.bins, CA)
#********************************************************************************************
#*  	Preprocess input bins as passed by the user or use default                          *
#********************************************************************************************
{
	if (class(CA$x) == "numeric") 
		{bins <- floor(sqrt(length(CA$x)))}								#If a vector - take sqrt of the length	
		else
		{bins <- floor(sqrt(nrow(CA$x)))} 					#Set the default number of bins 
		
	if (!is.na(input.bins))								#If User entered number of Bins
		if (suppressWarnings(!is.na(as.numeric(input.bins))))	#check if the User entered a numeric number of bins
			{	bins = input.bins	}					#Valid Input from the User - Use it		
			else
			{	cat('\nInvalid number of bins entered - Using default =nrow(x)) : ',bins,'\n')}
		
	CA$bins = bins										#Post it to Common Area
	return(CA)
}