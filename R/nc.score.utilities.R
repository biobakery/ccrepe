nc.score.helper <- function(data,CA) {
	#****************************************************************************************
	#*   Core calculations                                                                  *
	#****************************************************************************************
	mode(data) <- "numeric"
	CA$nc.score.matrix <- matrix(nrow=ncol(data),ncol=ncol(data))		#Build results area
	n <- length(unique(c(data)))
	adj <- ((1.5)*n*(n-1)/(n^2-n+1))
	for(i in 1:(ncol(data))) {    
		for(j in (i):ncol(data)) {   		
			ijsum <- 0
			cosum <- 0
			cesum <- 0
			for (m in 1:(nrow(data))) {   
				for (n in (m):nrow(data)) {  	  
					mx <- (c(data[m,i], data[m,j], data[n,i], data[n,i]))
         			if (length(unique(mx)) >= 2) {
						if (((mx[1]<mx[3])&(mx[2]<mx[4])) | ((mx[1]>mx[3])&(mx[2]>mx[4]))) {
							cosum <- cosum + 1 }
						if (((mx[1]>mx[3])&(mx[3]<mx[4])&(mx[4]>mx[2])&(mx[2]<mx[1])) | ((mx[1]<mx[3])&(mx[3]>mx[4])&(mx[4]<mx[2])&(mx[2]>mx[1]))) {
							cesum <- cesum + 1 }  
					}
				}
			}
		cesum_adj <- cesum * adj
		ijsum <- (cosum - cesum_adj)
		CA$nc.score.matrix[i,j] <- ijsum				#Post it to the matrix
		CA$nc.score.matrix[j,i] <- ijsum				#Post it to the matrix
		}
	}
  return(CA$nc.score.matrix)							#Return the calculated results matrix
}
nc.score.renormalize <-
function(
#*************************************************************************************
#*  	nc.score.renormalize                                                         *
#*************************************************************************************
   	x,						#First discretized entity 
	y,						#Second discretized entity
	nc.score.result)					  
{
	NS	<- nrow(as.matrix(x))										#Number of Samples
	NB 	<-   length(unique(c(x,y)))									#Number of bins
	renormalization.factor <-  choose(NS, 2) - (NS %% NB) * choose((floor(NS/NB) + 1), 2) - (NB - NS %% NB) * choose(floor(NS/NB), 2)
	if  (renormalization.factor == 0) 								#So that we dont get NaNs
		{renormalization.factor =  1} 		
	nc.score.result <- nc.score.result  / renormalization.factor	#Normalize
	return(nc.score.result)
}
nc.score.vectors.helper <-
function(
#*************************************************************************************
#*  	nc.score.helper                                                              *
#*************************************************************************************
   	x,						#First discretized  input 
	y)						#Second discretized input 
{
	ijsum <- 0				#Reset ijsum
	cosum <- 0				#Reset cosum
	cesum <- 0				#Reset cesum
	n <- length(unique(c(x)))	#<------ Need to QA and verify this!!! GW
	adj <- ((1.5)*n*(n-1)/(n^2-n+1))
	for (i in 1:(length(x)-1)) 		#Loop on the entries of x
		{
		for (j in (i+1):(length(y))) #Loop on the entries of y	
			{
				mx <- (c(x[i], y[i], x[j], y[j]))	# We have x and y - this is how mx looks in this conf. <---Need to verify !!!!! GW
				if (length(unique(mx)) >= 2) 							#These are the core calculations 
					{
					if (((mx[1]<mx[3])&(mx[2]<mx[4])) | ((mx[1]>mx[3])&(mx[2]>mx[4]))) 
						{ cosum <- cosum + 1 }
					if (((mx[1]>mx[3])&(mx[3]<mx[4])&(mx[4]>mx[2])&(mx[2]<mx[1])) | ((mx[1]<mx[3])&(mx[3]>mx[4])&(mx[4]<mx[2])&(mx[2]>mx[1]))) 
						{ cesum <- cesum + 1 }  
					}

			}
		}
	cesum_adj <- cesum * adj			
	ijsum <- (cosum - cesum_adj)
	return(ijsum)				
}
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
qc_filter <-
function(data,CA) {
#********************************************************************************************
#*  	quality control: select species with >= 1E-4 relative abundance in >= 10% of samples*
#********************************************************************************************
	tmp <- {}
	names <- {}

	for (i in 1:ncol(data)) 
	{
		if (length(which(data[,i] >= CA$min.abundance)) >= CA$min.samples*nrow(data))
		{
			tmp <- cbind(tmp, data[,i])
			names <- c(names, colnames(data)[i])
		}
	}
  colnames(tmp) <- names					#Post the column names
  rownames(tmp) <- rownames(data)			#Post the data row names
  return(tmp)
}