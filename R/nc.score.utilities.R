nc.score.helper <- function(data,CA) {
	#****************************************************************************************
	#*   Core calculations                                                                  *
	#****************************************************************************************
	mode(data) <- "numeric"
	CA$nc.score.matrix <- matrix(nrow=ncol(data),ncol=ncol(data))		#Build results area
	
	#***********************************************************
	#*  Vectorizing the calculations                           *
	#*  Calculating  looking at the matrix as a collection     *
	#*   of n_columns column vectors and performing the        *
	#*   calculations for the appropriate columns              *
	#***********************************************************
	

	for( i in seq_len(ncol(data)) ) 
	{  
		for( j in (i):(ncol(data)) ) 
		{
                        x <- data[,i]
			y <- data[,j]
                        n.bins <- length(unique(c(x,y)))
                        adj <- ((1.5)*n.bins*(n.bins-1)/(n.bins^2-n.bins+1))
			ri0 <- seq_len(length(x))
			rj0 <- Map(seq, ri0, length(y))
			ri <- rep(ri0, sapply(rj0, length))
			rj <- unlist(rj0)
			cosum <- sum(((x[ri]<x[rj])&(y[ri]<y[rj])) | ((x[ri]>x[rj])&(y[ri]>y[rj])))
			cesum <- sum(((x[ri]>x[rj])&(x[rj]<y[rj])& (y[rj]>y[ri])&(y[ri]<x[ri])) | ((x[ri]<x[rj])&(x[rj]>y[rj])& (y[rj]<y[ri])&(y[ri]>x[ri])))
			cesum_adj <- cesum * adj			
			ijsum <- (cosum - cesum_adj)
                        score <- nc.score.renormalize(x,y,ijsum)                        #Normalize the result
			CA$nc.score.matrix[i,j] <- score				#Post it to the matrix
			CA$nc.score.matrix[j,i] <- score				#Post it to the matrix
		}
	}
	return(CA$nc.score.matrix)
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
	n <- length(unique(c(x,y)))
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
        score <- nc.score.renormalize(x,y,ijsum)
	return(score)				
}
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
