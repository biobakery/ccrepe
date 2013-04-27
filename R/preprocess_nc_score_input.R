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


	if (nrow(x) == 1 & is.null(nrow(y)))					#If user entered x= vector and did not enter y
		{
		stop('Selected nc.score(frmOne=X) - X has to contain more than one row for meaningful results') 
		}
	if (is.null(x))									#If user entered x= vector and did not enter y
		{
		stop('Selected nc.score(frmOne=X) - X is null') 
		}
	
	CA <- list()
	Output <-list()
	x1 <- na.omit(x)								#Remove NAs from input rows
	CA <-list(data1 = t(x1),					    #Use the convention:Rows=Samples, cols=bugs)
			Output = Output,
			min.abundance = min.abundance,
			min.samples = min.samples
			)

	CA$YEntered = FALSE								#Default is that user did not enter Y
 
	if  (!is.null(nrow(y)))	
	   {								            #Check if user wants to check two vectors
		CA$YEntered = TRUE							#Set up the flag
		y1 <- na.omit(y)							#Remove NAs from input
		CA$data2 <- t(y1)							#Post Y1 into t(data2)
		if (!ncol(CA$data1) == 1   || !ncol(CA$data2) == 1 || !nrow(CA$data1) == nrow(CA$data2))	#Only allow vectors in this case
			{stop('For nc.score (x,y) x and y must be vectors of the same number of columns')}
		datamrg <- cbind(CA$data1,CA$data2)			#If they passed 2 parms
		CA$data1 <- datamrg							#We concatenate the vectors
		}
		
	

	CA$data1_trimmed <- qc_filter(CA$data1, CA)		#Trim the data
	
	CA$adBins <- sqrt(NROW(CA$data1)) 		#Set the default number of bins
	if (!is.na(adBins))								#If User entered number of Bins
		if (suppressWarnings(!is.na(as.numeric(adBins))))	#check if the User entered a numeric number of bins
			{	CA$adBins = adBins	}					#Valid Input from the User - Use it		
			else
			{	cat('\nInvalid number of bins entered - Using default =SQRT(nrow(data)) : ',CA$adBins,'\n')}

		
	if (!suppressWarnings(!is.na(as.numeric(CA$min.abundance))))	#check if the User entered a valid threshold1
		CA$min.abundance = 0.0001				#If it is not valid - force default
	if (!suppressWarnings(!is.na(as.numeric(CA$min.samples))))	#check if the User entered a valid threshold1
		CA$min.samples = 0.1				#If it is not valid - force default


	if (CA$adBins < 4)						#4 is the minimum number of bins
			{CA$adBins = 4}
	CA$data1_trimmed_cat <- discretize(CA$data1_trimmed,nbins=CA$adBins)	#Quantize relative abundance  data to ordinal categorical values 
	rownames(CA$data1_trimmed_cat) <- rownames(CA$data1_trimmed) #Discretize dropped them....			
	CA$data1_trimmed_cat_matrix <- as.matrix(CA$data1_trimmed_cat) #Convert the dataframe into matrix to process in NCSscore
	colnames(CA$data1_trimmed_cat_matrix) <- NULL
	return(CA)
}
