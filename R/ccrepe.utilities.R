calculate_q_values <-
function(CA, Gamma)
#**********************************************************************
#	Calculate Q values    				                              *
#**********************************************************************
{
	m = length(CA$p.values)					#m is the total number of p-values
	ln_m = log(m)									#log m
	ln_m_Plus_Gamma = ln_m + Gamma					
	SortedVector = sort(CA$p.values,index.return = TRUE)	#Sort the entire PValues matrix into a vector
	KVector = seq(1,m)						#A vector with the total number of entries in the PValues matrix
	QValues = SortedVector$x*m*ln_m_Plus_Gamma/KVector		#Calculate a vector containing the Q values
	QValuesArranged = rep(-1,m)
	QValuesArranged[SortedVector$ix] = QValues
	CA$q.values<-matrix(QValuesArranged, nrow= nrow(CA$p.values), byrow=TRUE)
	return(CA)
}




ccrepe_norm <-
function(data){
#*************************************************************************************
#* A function to normalize a data matrix across the rows                             *
#*************************************************************************************
	data.normalized <- data/apply(data,1,sum)		#Normalize
	data.normalized[is.na(data.normalized )] <- 0 #Replace NA's with Zeros
  	return(data.normalized)
}





ccrepe_process_one_dataset <-
function(data,N.rand, CA){
#*************************************************************************************
#* 	Function to apply the new method to a dataset                                    *
#*  As a note, data needs to be a matrix with no missing values                      *
#*************************************************************************************

	n = ncol(data)					# Number of columns, starting at 1; this is also the number of bugs
	
	CA$p.values <-matrix(data=0,nrow=n,ncol=n)	#Build the empty PValues matrix
	CA$cor <-matrix(data=0,nrow=n,ncol=n)	#Build the empty correlation matrix
	
	nsubj = nrow(data)				# Number of subjects

	# Generating the output data matrix (I chose this form for simplicity; we do want to keep the output object George has made)
	data.cor = matrix(ncol=5,nrow=choose(n,2))
	colnames(data.cor) = c("bug1", "bug2", "cor", "p.value", "q.value")

	# The matrix of possible bootstrap rows (which when multiplied by data give a specific row); of the form with all 0s except for one 1 in each row
	possible.rows = diag(rep(1,nsubj)) 

	#Generating the bootstrap and permutation matrices
	bootstrap.matrices = list()    # The list of matrices of possible bootstrap rows (each matrix multiplies by the data to give a resampled dataset)
	permutation.matrices = list()  # The list of permutation matrices; each matrix will be have columns which are permutations of the row indices

	for(i in 1:N.rand){
		if (CA$verbose == TRUE)		#If output is verbose - print the number of iterations
			{
			if (i %% CA$iterations.gap == 0)		#print status
				{cat('Completed ',i,' iterations','\n')}
			}
   
		# Get the rows of the possible.rows matrix; these correspond to the rows which will be included in the resampled dataset
		boot.rowids = sample(seq(1,nsubj),nsubj,replace=TRUE)

		# Generate the bootstrap matrix by getting the appropriate rows from the possible bootstrap rows
		boot.matrix = possible.rows[boot.rowids,] 

		# Add the bootstrap matrix to the list                   
		bootstrap.matrices = lappend(bootstrap.matrices,boot.matrix) 

		perm.matrix = replicate(n,sample(seq(1,nsubj),nsubj,replace=FALSE))  # The matrix has each column be a permutation of the row indices
		permutation.matrices = lappend(permutation.matrices,perm.matrix)     # Add the new matrix to the list
	}

 

	# The bootstrapped data; resample the data using each bootstrap matrix
	boot.data = lapply(bootstrap.matrices,resample,data=data)
	if (!is.null(CA$method.args.method))						#If user provided method.args.method use - otherwise don't
		{boot.cor = lapply(boot.data, CA$method, method=CA$method.args.method)}
		else
		{boot.cor = lapply(boot.data, CA$method)}

	# Generating the permutation data; permute the data using each permutation matrix
	permutation.data = lapply(permutation.matrices,permute,data=data)
	# The correlation matrices for the bootstrapped data; calculate the correlation for each resampled dataset

	# The normalized permutation data; the permutation data needs to be renormalized, but not the bootstrap data
	permutation.norm = lapply(permutation.data,ccrepe_norm)

	# The correlation matrices of the permuted data; calculate the correlation for each permuted dataset

	
	if (!is.null(CA$method.args.method))						#If user provided method.args.method use - otherwise don't
		{permutation.cor = lapply(permutation.norm,CA$method, method=CA$method.args.method)}
		else
		{permutation.cor = lapply(permutation.norm,CA$method) }
	# Now, actually calculating the correlation p-values within the dataset
    n.c = 0	# Counter for the number of comparisons (to enter in the output matrix)
	
	for(i in 1:n){
		if((i+1)<=n){
			for(k in (i+1):n){
				# Get a vector of the (i,k)th element of each correlation matrix in the list of bootstrapped data; this is the bootstrap distribution
				bootstrap.dist = unlist(lapply(boot.cor,'[',i,k)) 
          
				# Get a vector the (i,k)th element of each correlation matrix in the list of permuted data; this is the permuted distribution
				permutation.dist = unlist(lapply(permutation.cor,'[',i,k))	#sets
				
				n.0_1 = sum(data[,i]==0)				#Number of zeros in column i
				n.0_2 = sum(data[,k]==0)				#Number of zeros in column k
				n.p   = nrow(data)						#Number of rows in data
				CalcThresholdForError = ((CA$errthresh)^(1/n.p))*n.p		#If there is not enough data 
				if (n.0_1 > CalcThresholdForError | n.0_2 > CalcThresholdForError)
					{	
						p.value=NA
						cor=NA
					} else
					{
					if (!is.null(CA$method.args.method))	#If User provided method.args.method - use it -else - ignore it
						{cor = CA$method(data[,i],data[,k], method=CA$method.args.method)}   # Calculate the  correlation between the bugs
						else
						{cor = CA$method(data[,i],data[,k] )}  # Calculate the  correlation between the bugs
					# The Z-test to get the p-value for this comparison; as in get.renorm.null.pval
					
					
					 	
					
					p.value = pnorm(mean(permutation.dist), 
                                        mean=mean(bootstrap.dist), 
                                        sd=sqrt((var(bootstrap.dist) + var(permutation.dist))*0.5))
					}
				CA$p.values[i,k] = p.value				#Post it in the p-values matrix  
				CA$p.values[k,i] = p.value				#Post it in the p-values matrix  
				CA$cor[i,k] = cor						#Post it in the cor matrix  
				CA$cor[k,i] = cor						#Post it in the cor matrix  				
				n.c = n.c + 1
				data.cor[n.c,] = c(i,k,cor,p.value,NA)
			}
		}
	}
	
	
	
	
	Gamma = 0.57721566490153286060651209008240243104215933593992  		#I need to find the R version of Gamma!
	CA <- calculate_q_values(CA, Gamma)					#Calculate the QValues
	for (indx in 1:nrow(data.cor))						#post the q-values
		{
			i = data.cor[indx,1]
			k = data.cor[indx,2]
			data.cor[indx,5] = CA$q.values[i,k]
		}
	CA$data.cor <- data.cor								# Post it in the common Area

	if    (!is.na(CA$outdist))													#If user requested to print the distributions
		{
		lapply(boot.data,printDF,CA=CA,DistributionType=" Boot")
		lapply(permutation.norm,printDF,CA=CA,DistributionType=" Permutation")
		}
	if   (!is.na(CA$outdist))										#If user requested to print the distributions
		{
		close(CA$outdistFile)										#Close outdist file	
		CA$outdistFile <- NULL										#And remove it from the common area
		}
		
	#********************************************************************
	#*  Final Edits before exiting                                      *
	#********************************************************************
	diag(CA$p.values) <- NA											#Set diagonal of p.values to NA 
	colnames(CA$p.values)<-colnames(CA$data1)						#Set the names of the columns in the p.values matrix
	rownames(CA$p.values)<-colnames(CA$data1)						#Set the names of the roes in the p.values matrix
	colnames(CA$q.values)<-colnames(CA$data1)						#Set the names of the columns in the q.values matrix
	rownames(CA$q.values)<-colnames(CA$data1)						#Set the names of the roes in the q.values matrix
	colnames(CA$cor)<-colnames(CA$data1)							#Set the names of the columns in the q.values matrix
	rownames(CA$cor)<-colnames(CA$data1)							#Set the names of the roes in the q.values matrix
	return(CA)														# Return the output matrix
}





ccrepe_process_two_datasets <- function(data1.norm,data2.norm,N.rand, CA)
#*************************************************************************************
#* 	ccrepe function for two datasets                                                 *
#*************************************************************************************
{

	# Get number of bugs, subjects for each dataset
	n1 = ncol(data1.norm)
	n2 = ncol(data2.norm)

	data = merge_two_matrices(data1.norm,data2.norm)
	data1 <- data[,1:n1]
	data2 <- data[,(n1+1):(n1+n2)]

	
	
	nsubj = nrow(data)
	nsubj1 = nrow(data1)
	nsubj2 = nrow(data2)

	# Generating the output data matrix (I chose this form for simplicity; we do want to keep the output object George has made)
    bug1 <- rep(NA,n1*n2)
	bug2 <- rep(NA,n1*n2)
	cor.meas   <- rep(NA,n1*n2)
	p.values   <- rep(NA,n1*n2)

	# The matrix of possible bootstrap rows (which when multiplied by data give a specific row); of the form with all 0s except for one 1 in each row
	possible.rows = diag(rep(1,nsubj))

	#Generating the bootstrap and permutation matrices
	bootstrap.matrices   = list()  # The list of matrices of possible bootstrap rows (each matrix multiplies by the data to give a resampled dataset)
	permutation.matrices1 = list()  # The list of permutation matrices; each matrix will be have columns which are permutations of the row indices
	permutation.matrices2 = list()  # The list of permutation matrices; each matrix will be have columns which are permutations of the row indices

	for(i in 1:N.rand){

		# Get the rows of the two possible.rows matrices; these correspond to the rows which will be included in the resampled datasets
		boot.rowids        = sample(seq(1,nsubj),nsubj,replace=TRUE)

		# Generate the bootstrap matrices by getting the appropriate rows from the possible bootstrap rows
		boot.matrix        = possible.rows[boot.rowids,]

		# Add the bootstrap matrices to the appropriate lists  
		bootstrap.matrices = lappend(bootstrap.matrices,boot.matrix) # Add the bootstrap matrix to the list

		perm.matrix1 = replicate(n1,sample(seq(1,nsubj1),nsubj1,replace=FALSE)) # The matrix has each column be a permutation of the row indices
		permutation.matrices1 = lappend(permutation.matrices1,perm.matrix1)     # Add the new matrix to the list
		perm.matrix2 = replicate(n2,sample(seq(1,nsubj2),nsubj2,replace=FALSE)) # The matrix has each column be a permutation of the row indices
		permutation.matrices2 = lappend(permutation.matrices2,perm.matrix2)     # Add the new matrix to the list
	}


	# The bootstrapped data; resample the data using each bootstrap matrix
	boot.data = lapply(bootstrap.matrices,resample,data=data)

	# Generating the permutation data; permute the data using each permutation matrix
	permutation.data1 = lapply(permutation.matrices1,permute,data=data1)
	permutation.data2 = lapply(permutation.matrices2,permute,data=data2)

	# The normalized permutation data; the permutation data needs to be renormalized, but not the bootstrap data
	permutation.norm1 = lapply(permutation.data1,ccrepe_norm)
	permutation.norm2 = lapply(permutation.data2,ccrepe_norm)


	# The correlation matrices of the permuted and bootstrapped data; see extractCor function for more details
	# mapply is a function that applies over two lists, applying extractCor to the first element of each, then the
	# second element of each, etc.
 

		
	permutation.cor = mapply(extractCor,
		            mat1=permutation.norm1,
		            mat2=permutation.norm2,
		            startrow=1,
		            endrow=n1,
		            startcol=(n1+1),
		            endcol=(n1+n2), 				 
					use = "complete.obs",
		            method="spearman",
		            SIMPLIFY=FALSE)
		
			

	if (!is.null(CA$method.args.method))	#If User provided method.args.method - use it -else - ignore it
		{boot.cor = lapply(boot.data,CA$method,use="complete.obs",method=CA$method.args.method)}
		else
		{boot.cor = lapply(boot.data,CA$method,use="complete.obs")}

	



	# Now calculating the correlations and p-values between the two datasets
    n.c = 0	# Counter for the number of comparisons (to enter in the output matrix)
	
	
	CA$p.values <-matrix(data=0,nrow=n1,ncol=n2)	#Build the empty PValues matrix
	CA$cor <-matrix(data=0,nrow=n1,ncol=n2)	#Build the empty correlation matrix
	
	
	for(i in 1:n1){
		for(k in 1:n2){
			# Get a vector of the (i,k)th element of each correlation matrix in the list of bootstrapped data; this is the bootstrap distribution
			bootstrap.dist = unlist(lapply(boot.cor,'[',i,n1+k))
			bootstrap.dist[is.na(bootstrap.dist)] <- 0				#If there is an NA in bootstrap.dist - replace with 0 (Needs review)

			
			
			# Get a vector the (i,k)th element of each correlation matrix in the list of permuted data; this is the permuted distribution
			permutation.dist = unlist(lapply(permutation.cor,'[',i,k))

			n.c = n.c + 1
 			
			if (!is.null(CA$method.args.method))	#If User provided method.args.method - use it -else - ignore it
				{cor.meas[n.c] = CA$method(data[,i],data[,n1+k],use="complete.obs",method=CA$method.args.method)}
				else
				{cor.meas[n.c] = CA$method(data[,i],data[,n1+k],use="complete.obs")}
		
			
			
			# The Z-test to get the p-value for this comparison; as in get.renorm.null.pval
			p.value = pnorm(mean(permutation.dist), 
                                  mean=mean(bootstrap.dist), 
                                  sd=sqrt((var(bootstrap.dist) + var(permutation.dist))*0.5))

								  

								  
			CA$p.values[i,k] = p.value				#Post it in the p-values matrix  
			CA$cor[i,k] = cor.meas[n.c]				#Post it in the cor matrix  
							  
			p.values[n.c] = p.value
			bug1[n.c] = colnames(data1)[i]
			bug2[n.c] = colnames(data2)[k]
			

		}
	}
	#####data.cor <- data.frame(bug1,bug2,cor.meas,p.values)
	####return(data.cor)			# Return the output matrix
	

	CA$data.cor <- data.frame(bug1,bug2,cor.meas,p.values)
	rownames(CA$p.values) <- colnames(CA$data1.norm)		#Post the column names
	colnames(CA$p.values) <- colnames(CA$data2.norm)		#Post the column names
	rownames(CA$cor) <- colnames(CA$data1.norm)		#Post the column names
	colnames(CA$cor) <- colnames(CA$data2.norm)		#Post the column names
	CA$data1.norm <- NULL  # Not needed anymore
	CA$data2.norm <- NULL  # Not needed anymore
	CA$data1  <- NULL  # Not needed anymore
	CA$data2  <- NULL  # Not needed anymore
	CA$OneDataset <- NULL #Not needed anymore
	CA$subset.cols.1 <-NULL #Not needed anymore
	CA$subset.cols.2 <-NULL #Not needed anymore
	CA$method.args.method<- NULL #Not needed anymore
	CA$method <- NULL #Not needed anymore
	CA$verbose <- NULL #Not needed anymore
	CA$outdist <- NULL #Not needed
	CA$min.subj <- NULL #Not needed
	CA$errthresh <- NULL #Not needed
	CA$iterations.gap <-NULL #Not needed
	return(CA)			# Return the output matrix
}






ccrepe.calculations  <-
#*************************************************************************************
#*  	ccrepe calculations                                                          *
#*************************************************************************************
function(
 				x, 							#Data Frame  1 - For ccrepe and nc.score
				y,							#Data Frame  2 - For ccrepe and nc.score 
				method,						#Parameters,
				method.args.method,			#correlation method
				min.subj,					#Minimum rows in "data" frame to proceed (If not - run stops) - For ccrepe
				iterations,					#Reboot iterations - For ccrepe
				subset.cols.1,				#Subset of cols from cav1 to iterate on (c(0)== ALL) - For ccrepe
				subset.cols.2,				#Subset of cols from cav2 to iterate on (c(0)== ALL) - For ccrepe
				errthresh, 					#Threshold error if there is enough data to calculate cor an pval for certain i and k - For ccrepe
				verbose,					#Request for verbose output?
				iterations.gap,				#If output is verbose - after how many iterations issue a status message?
				distributions 				#Output Distribution file - For ccrepe
		)
{
	Output <-list()
	CA <-list(	data1=x,
			data2=y,
			min.subj=min.subj,
			iterations=iterations,
			errthresh=errthresh,
			method=method,
			method.args.method=method.args.method,
			subset.cols.1=subset.cols.1,
			subset.cols.2=subset.cols.2,
			verbose=verbose,
			iterations.gap=iterations.gap,
			outdist=distributions
			)

	
	CA <- DecodeInputCAArguments(CA)							#Decode Input Parameters

	if (CA$OneDataset == TRUE)
		{
		mydata.norm = preprocess_data(CA$data1,CA$subset.cols.1,CA)						#Preprocess the data 
		CA  = ccrepe_process_one_dataset(mydata.norm,CA$iterations, CA)  				#Invoke the new function 
		CA$data1 <- NULL										#Remove it from the output
		CA$data2 <- NULL										#Remove it from the output
		CA$subset.cols.2<-NULL									#Remove it from output
		CA$outdist <- NULL										#Remove it from output
		CA$OneDataset <- NULL									#Remove it from output
		CA$method <- NULL									#Remove it from output
		}
	else
		{
		CA$data1.norm = preprocess_data(CA$data1,CA$subset.cols.1,CA)					#Preprocess data1 
		CA$data2.norm = preprocess_data(CA$data2,CA$subset.cols.2,CA)					#Preprocess data2
		CA = ccrepe_process_two_datasets  (CA$data1.norm ,CA$data2.norm ,CA$iterations, CA)
		}
	return ( CA )
}






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
 
	if (identical(CA$method,cor) == TRUE)					#Check the methods argument only for cor
	{
	MethodsTable <- c('spearman','kendall','pearson')
	if   (is.null(CA$method.args.method) ||   CA$method.args.method  %in% MethodsTable == FALSE)
		{
		CA$method.args.method <- 'spearman'
		}
	}
	CA$data1 <- na.omit(CA$data1)						#Remove NA's
  
	if    (!is.na(CA$outdist))							#If the user passed a file - open it
		{
		CA$outdistFile = file(CA$outdist,open='at')		#Open outdist file
		}

	if ( CA$verbose !=  TRUE  & 	CA$verbose !=  FALSE )	#Verbose Flag has to be either true or false, otherwise - we set it to false
		{
		CA$verbose =  FALSE									#False - is the default
		}
	if  ( is.na(suppressWarnings(as.integer(CA$iterations.gap)))) 	#Check the iterations gap (Number of iterations after which to print status if verbose
		{ CA$iterations.gap = 100}						#If not valid - use 100
	return(CA)			 				#Return list of decoded input parameters
}




extractCor <-
function(mat1,mat2,startrow,endrow,startcol,endcol,method, ...)
#******************************************************************************************
# A function to calculate the correlation of the two matrices by merging them,            *
#     calculating the correlation of the merged matrix, and extracting the appropriate    *
#     submatrix to obtain the correlation of interest                                     *
# startrow, endrow, startcol, endcol give the indeces of the submatrix to extract         *
# The method refers to the correlation method; this function will need to be generalized  *
#******************************************************************************************
{
	mat <- merge_two_matrices(mat1,mat2)	             # Merge the two matrices
	mat_C <- cor(mat,method=method)                      # Calculate the correlation for the merged matrix
	sub_mat_C <- mat_C[startrow:endrow, startcol:endcol] # Extract the appropriate submatrix
	return(sub_mat_C)
}




merge_two_matrices <-
function(mat1,mat2)
#*************************************************************************************
# A function to merge two matrices, which works when they're different dimensions    *
# merge the two matrices by the rows, merging only rows present in both              *
# (ensures no missing samples)                                                       *
# exclude the first column, which is the row names                                   *
#*************************************************************************************
{
	mat <- as.matrix(merge(mat1,mat2,by="row.names"))
	rownames(mat) = mat[,1]    # The first column has the row names
	mat = mat[,-1]             # Remove the first column
	class(mat) <- "numeric"    # The default output of merge is strings; convert back to numbers
	return(mat)
}




lappend <-
function(lst, obj) {
#*************************************************************************************
#*  Function to append obj to lst, where lst is a list                               *
#*************************************************************************************
    lst[[length(lst)+1]] <- obj
    return(lst)
}





permute <-
function(data,permute.id.matrix){
#*************************************************************************************
#* Function to permute the data using a matrix of indices                            *
#* data is the data which you want to permute                                        *
#* permute.id.matrix is a matrix where each column is some                           *
#* permutation of the row indices                                                    * 
#*************************************************************************************
	permute.data = data
	for(i in 1:ncol(data)){									# For each column
		permute.data[,i] = data[permute.id.matrix[,i],i] 	# Replace the original column with a permuted one
	}
	return(permute.data)
}






preprocess_data <-
function(X,SelectedSubset,CA)
#**********************************************************************
#	Preprocess input data 				                              *
#**********************************************************************
{	
		MyDataFrame<-na.omit(X	)									#Post the input data into a working data frame - take only subset requestd by user
		MyDataFrame1  = MyDataFrame[,SelectedSubset]				#Select only the columns the User requested (If he did not: subset1=all columns)
		mydata <- MyDataFrame1[rowSums(MyDataFrame1 != 0) != 0, ] 	#Remove rows that are all zero to prevent NaNs
		if(nrow(mydata) < CA$min.subj) 						#If not enough data, issue messages in files and stop the run 
			{
			ErrMsg = paste('Not enough data - found ',nrow(mydata),' rows of data - Less rows than  ',CA$min.rows, ' min.rows - Run Stopped')  #Error 
			stop(ErrMsg)
			}
		mydata[1] <- NULL   #The first column is subject id - not data
		ProcessedX = mydata/apply(mydata,1,sum)
	return(ProcessedX)
}




printDF <-
function(Input1, CA, DistributionType)	{
#*************************************************************************************
#*  Function to print a data frame                                                   *
#*************************************************************************************
	OutputLine = c('Distribution:',DistributionType  )
	cat(OutputLine,file=CA$outdistFile,append=T)		#Output header to outdist file
	cat('\n',file=CA$outdistFile,append=T)				#Line feed
	df <- data.frame(Input1)							#convert input to data frame
	pm <- data.matrix(df, rownames.force = NA)			#convert df to matrix
	write.table(pm,CA$outdistFile,quote=F,row.names=F,col.names=F) 	#Write to file
	cat('\n',file=CA$outdistFile,append=T)		   		#Line feed
	return (0)
	}
	
	
	
	
	
	
resample <-
function(data,resample.matrix){
#*************************************************************************************
#* Function to calculate the resampled data matrix; this                             *
#*    is used for getting only the bootstrapped data                                 *
#*    data is the data which you want to bootstrap                                   *
#*   resample.matrix is the matrix appropriate for getting                           *
#*    you the bootstrapped dataset: matrix is square, with                           *
#*    the same number of rows as data; each row has exactly one                      *
#*    1 in it                                                                        *
#*************************************************************************************
	return(resample.matrix%*%data)
}
