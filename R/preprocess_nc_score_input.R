preprocess_nc_score_input <-
function(
#********************************************************************************************
#*  	Preprocess input and build the common area                                          *
#********************************************************************************************
			x, 										#First Input
			y,										#Second Input
			adBins,
			min.abundance,
			min.samples)
{
	if (ncol(x) == 1 & is.null(ncol(y)))					#If user entered x= vector and did not enter y
		{
		stop('Selected nc.score(frmOne=X) - X has to contain more than one row for meaningful results') 
		}
	if (is.null(x))									#If user entered x= vector and did not enter y
		{
		stop('Selected nc.score(frmOne=X) - X is null') 
		}
	
	CA <- list()
	Output <-list()
	CA <-list(data1 = t(x),					    		#Use the convention:Rows=Samples, cols=bugs)
			Output = Output,
			min.abundance = min.abundance,
			min.samples = min.samples
			)

	CA$YEntered = FALSE								#Default is that user did not enter Y
 
	if  (!is.null(nrow(y)))	
	   {								            #Check if user wants to check two vectors
		CA$YEntered = TRUE							#Set up the flag
		CA$data2 <- t(y)							#Post Y1 into t(data2)
		
		if (!nrow(CA$data1) == 1   || !nrow(CA$data2) == 1 || !ncol(CA$data1) == ncol(CA$data2))	#Only allow vectors in this case
			{stop('For nc.score (x,y) x and y must be vectors of the same number of columns')}
		datamrg <- rbind(CA$data1,CA$data2)			#If they passed 2 parms
		CA$data1 <- datamrg							#We concatenate the vectors
		}
		
		
	data1.transposed <- t(CA$data1)					#Transpose it to remove NAs
	data1.transposed.clean <- na.omit(data1.transposed)  # Remove rows that have NAs
	CA$data1.clean <-t(data1.transposed.clean) 				#Remove NAs from input rows
	

	CA$data1_trimmed <- qc_filter(CA$data1.clean, CA)		#Trim the data
	
	CA$adBins <- sqrt(ncol(CA$data1.clean)) 		#Set the default number of bins
	if (!is.na(adBins))								#If User entered number of Bins
		if (suppressWarnings(!is.na(as.numeric(adBins))))	#check if the User entered a numeric number of bins
			{	CA$adBins = adBins	}					#Valid Input from the User - Use it		
			else
			{	cat('\nInvalid number of bins entered - Using default =nrow(data)) : ',CA$adBins,'\n')}

		
	if (!suppressWarnings(!is.na(as.numeric(CA$min.abundance))))	#check if the User entered a valid threshold1
		{CA$min.abundance = 0.0001}				#If it is not valid - force default
	if (!suppressWarnings(!is.na(as.numeric(CA$min.samples))))	#check if the User entered a valid threshold1
		{CA$min.samples = 0.1}				#If it is not valid - force default

 
	CA$data1_trimmed_cat <- discretize(CA$data1_trimmed,nbins=CA$adBins)	#Quantize relative abundance  data to ordinal categorical values 
	rownames(CA$data1_trimmed_cat) <- rownames(CA$data1_trimmed) #Discretize dropped them....			
	CA$data1_trimmed_cat_matrix <- as.matrix(CA$data1_trimmed_cat) #Convert the dataframe into matrix to process in NCSscore
	colnames(CA$data1_trimmed_cat_matrix) <- NULL
	return(CA)
}
