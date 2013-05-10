preprocess_nc_score_input <-
function(
#********************************************************************************************
#*  	Preprocess input and build the common area                                          *
#********************************************************************************************
			x, 										#First Input
			y,										#Second Input
			input.adBins,
			input.min.abundance,
			input.min.samples)
{
 
	if (ncol(x) == 1 & is.null(ncol(y)))					#If user entered x= vector and did not enter y
		{
		stop('Selected nc.score(frmOne=X) - X has to contain more than one row for meaningful results') 
		}
	if (is.null(x))									#If user entered x= vector and did not enter y
		{stop('Please select y - ncscore(x,y)') }
	CA <-list()												#Common area definition 	 
 

 
 	if (!suppressWarnings(!is.na(as.numeric(input.min.abundance))))	#check if the User entered a valid threshold1
		{min.abundance = 0.0001}				#If it is not valid - force default
	CA$min.abundance = min.abundance

	if (!suppressWarnings(!is.na(as.numeric(input.min.samples))))	#check if the User entered a valid threshold1
		{min.samples = 0.1}				#If it is not valid - force default
	CA$min.samples = min.samples							#Post it to common Area
 
 
  
	#*********************************************************
	#* Filter and then select only the rows that are common  *
	#*********************************************************
	x <-qc_filter(x,CA)										#filter by abundance throsholds
	x <-na.omit(x)											#remove NAs
	x.keys.df <- as.data.frame(unlist(rownames(x)))			#To match
	colnames(x.keys.df)<-c('the.key')						#Give the column  a name		
	x.keys.df$included_x <- TRUE							#add flag
	
	
	y <-qc_filter(y,CA)										#filter by abundance throsholds
	y <-na.omit(y)											#remove NAs										 
	y.keys.df <- as.data.frame(unlist(rownames(y)))			#Build data frame of keys to match
	colnames(y.keys.df)<-c('the.key')						#Give the column  a name		
 	y.keys.df$included_y <- TRUE							#Add flag

	x.y.common.keys <- merge(x.keys.df, y.keys.df)			#Common Keys
	x.y.common.keys$included_x <- NULL						#Remove indicator - not needed anymore
	x.y.common.keys$included_y <- NULL						#Remove indicator - not needed anymore
	
	x.df <- as.data.frame(x)								#Build data frame from x
	x.df$the.key<-rownames(x.df)	 						#and build the key for merge
	x.common  <- merge(x.df, x.y.common.keys)				#Common Keys
	rownames(x.common) <- x.common$the.key 			
	x.common$the.key <- NULL								#Remove the key, not needed anymore
	
	y.df <- as.data.frame(y)								#Build data frame from x
	y.df$the.key<-rownames(y.df)	 						#and build the key for merge
	y.common  <- merge(y.df, x.y.common.keys)				#Common Keys
	rownames(y.common) <- y.common$the.key 					
	y.common$the.key <- NULL								#Remove the key, not needed anymore
	
	CA$x <- as.matrix(x.common)								#Post it to common Area
	CA$y <- as.matrix(y.common)								#Post it to common Area
	cat('\nAfter filtering there are ',nrow(CA$x),'common samples with ',ncol(CA$x),
		'variables (bugs) from the first input and ',ncol(CA$y),'variables (bugs) from the second input\n')
	CA$x.discretized <- matrix(nrow=nrow(CA$x),ncol=ncol(CA$x))				#Build x Discretized empty Matrix			
	CA$y.discretized <- matrix(nrow=nrow(CA$y),ncol=ncol(CA$y))				#Build y Discretized empty Matrix	
	#****************************************************************************************
	#*   Process Bins  passed by the user                                                   *
	#****************************************************************************************
	adBins <- floor(sqrt(nrow(CA$x))) 							#Set the default number of bins
	if (!is.na(input.adBins))								#If User entered number of Bins
		if (suppressWarnings(!is.na(as.numeric(input.adBins))))	#check if the User entered a numeric number of bins
			{	adBins = input.adBins	}					#Valid Input from the User - Use it		
			else
			{	cat('\nInvalid number of bins entered - Using default =nrow(x)) : ',adBins,'\n')}
		
	CA$adbins = adBins										#Post it to Common Area



	return(CA)
}
