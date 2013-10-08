nc.score.helper <- function(data,CA) {
	#****************************************************************************************
	#*   Core calculations                                                                  *
	#****************************************************************************************
	mode(data) <- "numeric"
	CA$nc.score.matrix <- matrix(nrow=ncol(data),ncol=ncol(data))		#Build results area
	n <- CA$n.bins
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
	y,						#Second discretized input 
	CA)						#Common area
{
	ijsum <- 0				#Reset ijsum
	cosum <- 0				#Reset cosum
	cesum <- 0				#Reset cesum
	n <- CA$n.bins
	adj <- ((1.5)*n*(n-1)/(n^2-n+1))    
	#**************************************************************
	# Vectorized the calculations                                 *
	#**************************************************************
	i0 <- seq_len(length(x) - 1)
    j0 <- Map(seq, i0 + 1, length(y))
    i <- rep(i0, sapply(j0, length))
    j <- unlist(j0)
    cosum <- sum(((x[i]<x[j])&(y[i]<y[j])) | ((x[i]>x[j])&(y[i]>y[j])))
    cesum <- sum(((x[i]>x[j])&(x[j]<y[j])& (y[j]>y[i])&(y[i]<x[i])) |   ((x[i]<x[j])&(x[j]>y[j])& (y[j]<y[i])&(y[i]>x[i])))
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

	for (i in 1:ncol(data)) 
	{
		if (length(which(data[,i] >= CA$min.abundance)) >= CA$min.samples*nrow(data))
		{
			tmp <- cbind(tmp, data[,i])
			names <- c(names, colnames(data)[i])
		} else 
		{
			warning_msg <- paste(colnames(data)[i],"doesn't pass quality control: try adjusing min.samples or min.abundance")
			warning(warning_msg)
		}
	}
  colnames(tmp) <- names					#Post the column names
  rownames(tmp) <- rownames(data)			#Post the data row names
  return(tmp)
}
