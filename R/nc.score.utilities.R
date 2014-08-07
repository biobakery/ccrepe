preprocess_nc_score_input <-
function(
#********************************************************************************************
#*  	Pre process input and build the common area                                         *
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
	CA$x <- x												#Post x in the common area
	CA <-qc_filter(x,CA)										#filter by abundance thresholds
	x <-na.omit(CA$x)											#remove NAs
	if (is.null(rownames(x))) {rownames(x)<-seq(1:nrow(x))} #If there are no row names - plug them in
	CA$x <-  x 												#Post it to common Area
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
	if (length(input.bins) > 1)					#If the user passed a set of bins,  sort them and use them
		{
			CA$bins <- sort(input.bins)
			CA$n.bins <- length(intersect(which(CA$bins<1),which(CA$bins>0)))
			if (! 0 %in% CA$bins )				#Check if 0 is in the list of the bins the User passed - not there - add it and sort 
				{ 
				 CA$bins<-c(CA$bins,0)			#Add 0
				 CA$bins <- sort(CA$bins)		#And sort it
				 CA$n.bins <- length(intersect(which(CA$bins<1),which(CA$bins>0)))  # The number of bins, since 1 and 0 represent the ends
				}
			return(CA)
		}
		
	if (class(CA$x) == "numeric") 
		{
		bins <- floor(sqrt(length(CA$x)))					#If a vector - take sqrt of the length	
		} else { 
		bins <- floor(sqrt(nrow(CA$x))) 					#Set the default number of bins 
		}
		
	if (!is.na(input.bins))								#If User entered number of Bins
		if (suppressWarnings(!is.na(as.numeric(input.bins))))	#check if the User entered a numeric number of bins
			{	
			bins = input.bins					#Valid Input from the User - Use it		
			} else 
			{
			warning('Invalid number of bins entered - Using default bins = sqrt(Number of rows) ')
			}


			
	CA$bins = bins										#Post it to Common Area
	CA$n.bins = bins									#Claculate the number of bins
	return(CA)
}
qc_filter <-
function(data,CA) {
#********************************************************************************************
#*  	quality control: select species with >= 1E-4 relative abundance in >= 10% of samples*
#********************************************************************************************
	tmp <- {}
	names <- {}
	CA$names.of.cols.failing.qc <- {} 			#Names of the cols failing QC
	CA$original.column.names <- colnames(data)	#These are the original column names
	
	CA$input.total.cols <- ncol(data)			#number of cols in the input
	CA$columns.not.passing.qc = vector()		#define a vector to contain the seq number of cols that did not pass qc

	

	for (i in seq_len(ncol(data))) 
	{
		if (length(which(data[,i] >= CA$min.abundance)) >= CA$min.samples*nrow(data))
		{
			tmp <- cbind(tmp, data[,i])
			names <- c(names, colnames(data)[i])
		} else 
		{
			CA$names.of.cols.failing.qc<-c(CA$names.of.cols.failing.qc, colnames(data)[i])
			CA$columns.not.passing.qc <- c(CA$columns.not.passing.qc,i)		#Add the number of the col that did not pass
			warning_msg <- paste("Column ",i,":  ",colnames(data)[i]," - doesn't pass quality control: try adjusting min.samples or min.abundance")
			warning(warning_msg)
		}
	}
  
  colnames(tmp) <- names					#Post the column names
  rownames(tmp) <- rownames(data)			#Post the data row names
  CA$x <- tmp								#Post in the common area
  return(CA)								#Return the common Area
}
