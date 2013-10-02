 calculate_q_values <- function(CA)
#**********************************************************************
#	Calculate Q values    				                              *
#**********************************************************************
{
	q.values.calc.matrix <- CA$p.values   #Set the default p.values  for one and two dataset cases
	if  (CA$OneDataset == TRUE)	
		{q.values.calc.matrix[lower.tri(q.values.calc.matrix)] <- NA}    #If One dataset,  use only the upper part of the matrix
	non.missing.p.values <- q.values.calc.matrix[which(!is.na(q.values.calc.matrix))]
	m = length(non.missing.p.values)                                #m is the total number of p-values
	ln_m = log(m)									#log m
	ln_m_Plus_Gamma = ln_m + CA$Gamma					
	SortedVector = sort(q.values.calc.matrix,index.return = TRUE)	#Sort the entire PValues matrix into a vector
	KVector = seq(1,m)						#A vector with the total number of entries in the PValues matrix
	QValues = SortedVector$x*m*ln_m_Plus_Gamma/KVector		#Calculate a vector containing the Q values
	QValuesArranged = rep(-1,m)
	QValuesArranged[SortedVector$ix] = QValues
	QValues = SortedVector$x*m*ln_m_Plus_Gamma/KVector		#Calculate a vector containing the Q values
	CA$q.values<-q.values.calc.matrix
	CA$q.values[which(!is.na(q.values.calc.matrix))] = QValuesArranged
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
	CA$z.stat <-matrix(data=0,nrow=n,ncol=n)	#Build the empty z.stat matrix
	
	CA$cor <-matrix(data=0,nrow=n,ncol=n)	#Build the empty correlation matrix
	
	nsubj = nrow(data)				# Number of subjects

	# Generating the output data matrix (I chose this form for simplicity; we do want to keep the output object George has made)
	
	data.cor = matrix(ncol=5,nrow=(n+1)*(n/2))	#Row Number increased from range(n,2) to n+1*(n/2) because we calc diagonal
	colnames(data.cor) = c("bug1", "bug2", "cor", "p.value", "q.value")

	# The matrix of possible bootstrap rows (which when multiplied by data give a specific row); of the form with all 0s except for one 1 in each row
	possible.rows = diag(rep(1,nsubj)) 

	#Generating the bootstrap and permutation matrices
	bootstrap.matrices = list()    # The list of matrices of possible bootstrap rows (each matrix multiplies by the data to give a resampled dataset)
	permutation.matrices = list()  # The list of permutation matrices; each matrix will be have columns which are permutations of the row indices
	
	

	if (length(colnames(CA$data1)) == 0)		#If no colum names - force them to be the column number
		{
		colnames(CA$data1)<-1:ncol(CA$data1) 
		}

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
	

	###  Using method.calculation  for the one dataset also
	boot.cor  = lapply(boot.data,method.calculation,nsubj,data,CA )		 ###Function to check is all cols are zeros and apply cor


	# Generating the permutation data; permute the data using each permutation matrix
	permutation.data = lapply(permutation.matrices,permute,data=data)
	# The correlation matrices for the bootstrapped data; calculate the correlation for each resampled dataset

	# The normalized permutation data; the permutation data needs to be renormalized, but not the bootstrap data
	permutation.norm = lapply(permutation.data,ccrepe_norm)

	# The correlation matrices of the permuted data; calculate the correlation for each permuted dataset
	
	
	CA$verbose.requested = FALSE			#If the User requested verbose output - turn it off temporarily
	if (CA$verbose == TRUE)
		{
			CA$verbose.requested <- TRUE		#Turn of the verbose flag save variable
			CA$verbose <- FALSE					#But pass non verbose to the CA$method (Could be nc.score or anything else)	
		}
	permutation.cor <- do.call(lapply,c(list(permutation.norm,CA$method), CA$method.args))  #Invoke the measuring function
	 
	# Now, actually calculating the correlation p-values within the dataset
    n.c = 0	# Counter for the number of comparisons (to enter in the output matrix)
	
 

	loop.range <- 1:n						#Establish looping range default
	if ( length(CA$subset.cols.1) > 1 )		#If the User entered a subset of columns
		{
		loop.range <- CA$subset.cols.1		#Use the subset of columns
		}


	for( index1 in 1:(length(loop.range)) )
	{
		i = loop.range[index1]
		{
			for(index2 in (index1):length(loop.range))
			{	
				k = loop.range[index2]
				# Get a vector of the (i,k)th element of each correlation matrix in the list of bootstrapped data; this is the bootstrap distribution
				
				bootstrap.dist = unlist(lapply(boot.cor,'[',i,k)) 
          
				# Get a vector the (i,k)th element of each correlation matrix in the list of permuted data; this is the permuted distribution
				permutation.dist = unlist(lapply(permutation.cor,'[',i,k))	#sets
				
				
				if    (!is.na(CA$outdist))						#If user requested to print the distributions
					{
					RC <- print.dist(bootstrap.dist,permutation.dist,CA,i,k)
					}
					
				n.0_1 = sum(data[,i]==0)				#Number of zeros in column i
				n.0_2 = sum(data[,k]==0)				#Number of zeros in column k
				n.p   = nrow(data)						#Number of rows in data

				CalcThresholdForError = ((CA$errthresh)^(1/n.p))*n.p		#If there is not enough data 
				if (n.0_1 > CalcThresholdForError | n.0_2 > CalcThresholdForError)
					{	
						p.value=NA
						cor=NA
                                                z.stat=NA
					} else
					{

					measure.parameter.list <- append(list(x=data[,i],y=data[,k]), CA$sim.score.parameters)  #build the method do.call parameter list
					cor <- do.call(CA$method,measure.parameter.list)	#Invoke the measuring function

					####################################################
					#  New p.value calculation                         #
					####################################################
					z.stat <- (mean(bootstrap.dist) - mean(permutation.dist))/sqrt(0.5*(var(permutation.dist)+var(bootstrap.dist)))
					p.value <- 2*pnorm(-abs(z.stat))					
					}
				CA$z.stat[i,k] = z.stat					#Post z.stat in output matrix	
				CA$z.stat[k,i] = z.stat					#Post z.stat in output matrix					
				CA$p.values[i,k] = p.value				#Post it in the p-values matrix  
				CA$p.values[k,i] = p.value				#Post it in the p-values matrix  
				CA$cor[i,k] = cor						#Post it in the cor matrix  
				CA$cor[k,i] = cor						#Post it in the cor matrix  				
				n.c = n.c + 1
				data.cor[n.c,] = c(i,k,cor,p.value,NA)
			}
		}
	}
	

	
	CA <- calculate_q_values(CA)						#Calculate the QValues
	CA$q.values[lower.tri(CA$q.values)] = CA$q.values[upper.tri(CA$q.values)]  #Making the q.values matrix symmetrical for the one dataset case
	
	for (indx in 1:nrow(data.cor))						#post the q-values
		{
			i = data.cor[indx,1]
			k = data.cor[indx,2]
			data.cor[indx,5] = CA$q.values[i,k]
		}
	CA$data.cor <- data.cor								# Post it in the common Area

	

		
	#********************************************************************
	#*  Final Edits before exiting                                      *
	#********************************************************************
	diag(CA$q.values) <- NA											#Set diagonal of p.values to NA 
	colnames(CA$p.values)<-colnames(CA$data1)						#Set the names of the columns in the p.values matrix
	rownames(CA$p.values)<-colnames(CA$data1)						#Set the names of the roes in the p.values matrix
	colnames(CA$q.values)<-colnames(CA$data1)						#Set the names of the columns in the q.values matrix
	rownames(CA$q.values)<-colnames(CA$data1)						#Set the names of the roes in the q.values matrix
	colnames(CA$cor)<-colnames(CA$data1)							#Set the names of the columns in the q.values matrix
	rownames(CA$cor)<-colnames(CA$data1)							#Set the names of the roes in the q.values matrix
	CA$sim.score <- CA$cor											#Rename cor to sim.score
	diag(CA$p.values) <- NA											#Set diagonal of p.values to NA 
	CA$cor <- NULL
	CA <- clean_common_area_after_processing(CA)	#Clean the Common Area before returning to the User
 
	if (length(CA$subset.cols.x) > 1)				#If used a subset - present only the subset
		{
		CA$p.values <- CA$p.values[CA$subset.cols.x,CA$subset.cols.x]   #Display only the subset of cols and rows
		CA$q.values <- CA$q.values[CA$subset.cols.x,CA$subset.cols.x]   #Display only the subset of cols and rows
		CA$sim.score <- CA$sim.score[CA$subset.cols.x,CA$subset.cols.x]   #Display only the subset of cols and rows
		}
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
	if(nrow(data) < CA$min.subj ) 	#If not enough data, issue messages in files and stop the run 
			{
			ErrMsg = paste('Not enough data - found ',nrow(data),' rows of data in the merged matrix - Less than  ',CA$min.subj, ' (=min.subj)  - Run Stopped')  #Error 
			stop(ErrMsg)
			}
	
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


	if (length(colnames(CA$data1.norm)) == 0)		#If no colum names - force them to be the column number
		{
		colnames(CA$data1.norm)<-1:ncol(CA$data1.norm) 
		}
	if (length(colnames(CA$data2.norm)) == 0)		#If no colum names - force them to be the column number
		{
		colnames(CA$data2.norm)<-1:ncol(CA$data2.norm) 
		}	
	
	
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
 
	CA$verbose.requested = FALSE			#If the User requested verbose output - turn it off temporarily
	if (CA$verbose == TRUE)
		{
			CA$verbose.requested <- TRUE		#Turn of the verbose flag save variable
			CA$verbose <- FALSE					#But pass non verbose to the CA$method (Could be nc.score or anything else)	
		}

	permutation.cor = mapply(extractCor,
		            mat1=permutation.norm1,
		            mat2=permutation.norm2,
		            startrow=1,
		            endrow=n1,
		            startcol=(n1+1),
		            endcol=(n1+n2), 				 
					MoreArgs = list(my.method=CA$method,
					method.args=CA$method.args,
					outdist=CA$outdist,
					outdistFile=CA$outdistFile),
		            SIMPLIFY=FALSE)
		
			

	#******************************************************************************** 
	# For each matrix, check if there is a column that is all zeros 				*
	# If such matrix is found,  try to reboot it at most 5 times until such matrix  *
	# is found that does not contail a column with all zeros                        *
	# The limit of 5 tries is now hard coded and needs to be re-evaluated           *
	# If after 5 times there is no success - we stop the run (Need to verify!!!!    *
	#********************************************************************************
	
	
	boot.cor  = lapply(boot.data,method.calculation,nsubj,data,CA )		 ###Function to check is all cols are zeros and apply cor

	
	# Now calculating the correlations and p-values between the two datasets
    n.c = 0	# Counter for the number of comparisons (to enter in the output matrix)
	
	
	CA$p.values <-matrix(data=0,nrow=n1,ncol=n2)	#Build the empty p.values matrix
	CA$z.stat  <-matrix(data=0,nrow=n1,ncol=n2)		#Build the empty z.stat matrix
	CA$cor <-matrix(data=0,nrow=n1,ncol=n2)	#Build the empty correlation matrix
	
	loop.range1 <- 1:n1						#Establish looping range default
	
	
	if ( length(CA$subset.cols.1)> 1 )		#If the User entered a subset of columns
		{
		loop.range1 <- CA$subset.cols.1		#Use the subset of columns
		}

	loop.range2 <- 1:n2						#Establish looping range default
	
	if ( length(CA$subset.cols.2) > 1 )		#If the User entered a subset of columns
		{
		loop.range2 <- CA$subset.cols.2		#Use the subset of columns
		}
	
	
	

	for(index1 in 1:(length(loop.range1)))
	{
		i = loop.range1[index1]
		for(index2 in 1:(length(loop.range2)))
		{
			k = loop.range2[index2]
			# Get a vector of the (i,k)th element of each correlation matrix in the list of bootstrapped data; this is the bootstrap distribution

			bootstrap.dist = unlist(lapply(boot.cor,'[',i,n1+k))

			bootstrap.dist[is.na(bootstrap.dist)] <- 0				#If there is an NA in bootstrap.dist - replace with 0 (Needs review)
			
			# Get a vector the (i,k)th element of each correlation matrix in the list of permuted data; this is the permuted distribution
			permutation.dist = unlist(lapply(permutation.cor,'[',i,k))
			
			if    (!is.na(CA$outdist))						#If user requested to print the distributions
					{
					RC <- print.dist(bootstrap.dist,permutation.dist,CA,i,k)
					}

			n.c = n.c + 1
		
					

			n.0_1 = sum(data[,i]==0)				#Number of zeros in column i
			n.0_2 = sum(data[,n1+k]==0)				#Number of zeros in column n1+k
			n.p   = nrow(data)						#Number of rows in data
			CalcThresholdForError1 = ((CA$errthresh)^(1/n.p))*n.p		#If there is not enough data on col1
			CalcThresholdForError2 = ((CA$errthresh)^(1/n.p))*n.p		#If there is not enough data on col1
		
	
			
			if (n.0_1 > CalcThresholdForError1 | n.0_2 > CalcThresholdForError2)
					{	
						p.value=NA
						cor=NA
                                                z.stat=NA
					} 
			else
					{
						measure.parameter.list <- append(list(x=data[,i],y=data[,n1+k]), CA$sim.score.parameters)  #build the method do.call parameter list
						cor.meas[n.c] <- do.call(CA$method,measure.parameter.list)	#Invoke the measuring function

						####################################################
						#  New p.value calculation                         #
						####################################################
						z.stat <- (mean(bootstrap.dist) - mean(permutation.dist))/sqrt(0.5*(var(permutation.dist)+var(bootstrap.dist)))
						p.value <- 2*pnorm(-abs(z.stat))					
					}

								  
			CA$p.values[i,k] = p.value				#Post it in the p.values matrix  
			CA$z.stat[i,k] = z.stat					#Post it in the z.stat matrix 
			CA$cor[i,k] = cor.meas[n.c]				#Post it in the cor matrix  
							  
			p.values[n.c] = p.value
			bug1[n.c] = colnames(data1)[i]
			bug2[n.c] = colnames(data2)[k]
			

		}
	}
	

	CA <- calculate_q_values(CA)						#Calculate the QValues

	CA$data.cor <- data.frame(bug1,bug2,cor.meas,p.values)
	
	for (indx in 1:nrow(CA$data.cor))						#post the q-values
		{
			i = CA$data.cor[indx,1]
			k = CA$data.cor[indx,2]
			CA$data.cor[indx,5] = CA$q.values[i,k]
		}
	colnames(CA$data.cor)[3] <- 'cor'						#Set up the name of the column
	colnames(CA$data.cor)[5] <- 'q.values'					#Set up the name of the column


	
	#********************************************************************
	#*  Final Edits before exiting                                      *
	#********************************************************************
	rownames(CA$p.values) <- colnames(CA$data1.norm)		#Post the column names
	colnames(CA$p.values) <- colnames(CA$data2.norm)		#Post the column names
	rownames(CA$cor) <- colnames(CA$data1.norm)		#Post the column names
	colnames(CA$cor) <- colnames(CA$data2.norm)		#Post the column names
	rownames(CA$q.values) <- colnames(CA$data1.norm)		#Post the column names
	colnames(CA$q.values) <- colnames(CA$data2.norm)		#Post the column names
	rownames(CA$z.stat) <- colnames(CA$data1.norm)		#Post the column names
	colnames(CA$z.stat) <- colnames(CA$data2.norm)		#Post the column names
	diag(CA$q.values) <- NA									#Set diagonal of q.values to NA 
	diag(CA$p.values) <- NA									#Set diagonal of p.values to NA 
	CA$sim.score <- CA$cor									#Rename cor to sim.score
	CA$cor <- NULL
	CA <- clean_common_area_after_processing(CA)	#Clean the Common Area before returning to the User
	total.rows.to.display = 1:nrow(CA$p.values)				#Number of rows to display
	total.cols.to.display = 1:ncol(CA$p.values)				#Number of cols to display

	if  (length(CA$subset.cols.x) > 1) #If User selected a subset
		{ 
			total.rows.to.display = CA$subset.cols.x
		}
	if  (length(CA$subset.cols.y) > 1) #If User selected a subset
		{ 
			total.cols.to.display = CA$subset.cols.y
		}
	CA$p.values <- CA$p.values[total.rows.to.display,total.cols.to.display]   #Display only the subset of cols and rows
	CA$q.values <- CA$q.values[total.rows.to.display,total.cols.to.display]   #Display only the subset of cols and rows
	CA$sim.score <- CA$sim.score[total.rows.to.display,total.cols.to.display]   #Display only the subset of cols and rows
	return(CA)			# Return the output matrix
}






ccrepe.calculations  <-
#*************************************************************************************
#*  	ccrepe calculations                                                          *
#*************************************************************************************
function(CA) 
{
	Output <-list()
 
	CA <- DecodeInputCAArguments(CA)							#Decode Input Parameters

	if (CA$OneDataset == TRUE)
		{
		mydata.norm = preprocess_data(CA$data1,CA)						#Preprocess the data 
		CA  = ccrepe_process_one_dataset(mydata.norm,CA$iterations, CA)  				#Invoke the new function 
		}
	else
		{
		CA$data1.norm = preprocess_data(CA$data1,CA)					#Preprocess data1 
		CA$data2.norm = preprocess_data(CA$data2,CA)					#Preprocess data2
		CA = ccrepe_process_two_datasets  (CA$data1.norm ,CA$data2.norm ,CA$iterations, CA)
		}
	return ( CA )
}






DecodeInputCAArguments <-
function(CA){
#*************************************************************************************
#*	Decode and validate input parameters                                             *
#*************************************************************************************
	CA$Gamma = 0.57721566490153286060651209008240243104215933593992  		#Euler's Gamma


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

 
 
	if  ( length(CA$subset.cols.1) == 0 && (is.null(CA$subset.cols.1)))				#If NULL- set to default
		{
		CA$subset.cols.1 = c(0)					#Set to default
		}
	if  ( length(CA$subset.cols.2) == 0 && (is.null(CA$subset.cols.2)))				#If NULL - set to default
		{
		CA$subset.cols.2 = c(0)					#Set to default
		}
	
  
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
		

	if  (identical(cor,CA$method) && length(CA$method.args) == 0)	#If the method is cor and the User did not pass any parms
		{
			CA$method.args = list(method='spearman',use='complete.obs')		#Set the default for cor
		}
	for (name in names(CA$method.args)) {					#Add the entries in method.args to the measuring parameter list			
 		CA$sim.score.parameters[[name]]<-CA$method.args[[name]]
		}	
		

	CA$retries.max.iterations =  -round(log2(CA$errthresh)) #This is the maximum number of iterations to try to reboot a matrix if in a col all values are 0
		
		
	return(CA)			 				#Return list of decoded input parameters
}




extractCor <-
function(mat1,mat2,startrow,endrow,startcol,endcol,my.method,method.args,outdist,outdistFile,  ...)
#******************************************************************************************
# A function to calculate the correlation of the two matrices by merging them,            *
#     calculating the correlation of the merged matrix, and extracting the appropriate    *
#     submatrix to obtain the correlation of interest                                     *
# startrow, endrow, startcol, endcol give the indeces of the submatrix to extract         *
# The method and the method parameters are passed via the list
#******************************************************************************************
{
  	mat <- merge_two_matrices(mat1,mat2)	            #Merge the two matrices
	measure.function.parm.list <- append(list(x=mat), method.args)	
	mat_C <-do.call(my.method,measure.function.parm.list)	#Invoke the measuring fnction
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
# <----------- Important  !!! --------------->                                       *
#* The matrices are merged by row name (Which defaults to row number)                *
#* so the assumption is that SAME ROWS pertain to SAME SUBJECTS                      *
#* Also not that first column is treated as data - not Subject.ID or any other       *
#* identifier                                                                        *
#*************************************************************************************
{
 	mat <- merge(mat1,mat2,by="row.names")
	rownames(mat) = mat[,1]    # The first column has the row names
	mat = mat[,-1]             # Remove the first column
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
function(X,CA)
#**********************************************************************
#	Preprocess input data 				                              *
#**********************************************************************
{	
		MyDataFrame<-na.omit(X	)									
		mydata <- MyDataFrame[rowSums(MyDataFrame != 0) != 0, ] 	#Remove rows that are all zero to prevent NaNs
		if(nrow(mydata) < CA$min.subj && CA$OneDataset == TRUE ) 	#If not enough data, issue messages in files and stop the run 
			{
			ErrMsg = paste('Not enough data - found ',nrow(mydata),' rows of data - Less than  ',CA$min.subj, ' (=min.subj) - Run Stopped')  #Error 
			stop(ErrMsg)
			}
		ProcessedX = mydata/apply(mydata,1,sum)
	return(ProcessedX)
}




	
print.dist <-
function(bootstrap.dist,permutation.dist,  CA,i, k)	{
#*************************************************************************************
#*  Function to print the distribution                                               *
#*************************************************************************************
	outdistFile = CA$outdistFile 						#Output file
	
	output.string0 <- paste(bootstrap.dist,sep="",collapse=',')		#Build the output string 
	if (CA$OneDataset == TRUE)							#The structure of the output is different for one dataset and two datasets
		{output.string <-paste("Boot,",colnames(CA$data1)[i],",",colnames(CA$data1)[k],",",output.string0,'\n', sep='',collapse=",")}
		else
		{output.string <-paste("Boot,",colnames(CA$data1.norm)[i],",",colnames(CA$data2.norm)[k],",",output.string0,'\n', sep='',collapse=",")}

	cat(output.string,file=CA$outdistFile,append=TRUE)
	
	output.string0 <- paste(permutation.dist,sep="",collapse=',')		#Build the output string 
	if (CA$OneDataset == TRUE)							#The structure of the output is different for one dataset and two datasets
		{output.string <-paste("Permutation,",colnames(CA$data1)[i],",",colnames(CA$data1)[k],",",output.string0,'\n', sep='',collapse=",")}
		else
		{output.string <-paste("Permutation,",colnames(CA$data1.norm)[i],",",colnames(CA$data2.norm)[k],",",output.string0,'\n', sep='',collapse=",")}
	cat(output.string,file=CA$outdistFile,append=TRUE)


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
	data.matrix <-as.matrix(data)
	return(resample.matrix%*%data.matrix)
}




clean_common_area_after_processing <-
function(CA){
#*************************************************************************************
#* 	Function to clean the Common Area After processing                               *
#*  Common Area is passed to the User as results                                     *
#*************************************************************************************
	if   (!is.na(CA$outdist))										#If user requested to print the distributions
		{
		close(CA$outdistFile)										#Close outdist file	
		CA$outdistFile <- NULL										#And remove it from the common area
		}
	if (CA$verbose.requested == TRUE)								#Check if the User requested verbose
		{															#We might have turned it off before calling CA$method	
		CA$verbose = TRUE	
		}
	CA$verbose.requested = NULL		
	CA$data1 <- NULL													
	CA$data1.norm <- NULL  											
	CA$data2 <- NULL													
	CA$data2.norm <- NULL  												
	CA$method <- NULL													
	CA$method.args<- NULL 										
	CA$OneDataset <- NULL													
	CA$outdist <- NULL	
	CA$Gamma <- NULL 
	CA$data.cor <- NULL 
	CA$retries.max.iterations <- NULL
	CA$subset.cols.x <-CA$subset.cols.1
	CA$subset.cols.y <-CA$subset.cols.2
	CA$subset.cols.1 <- NULL
	CA$subset.cols.2 <- NULL
	

	if  ( length(CA$subset.cols.x) == 1 &&  CA$subset.cols.x == c(0))			#If NA - set to default
		{
		CA$subset.cols.x <-NULL					#Set to default
		}
	if  ( length(CA$subset.cols.y) == 1 &&  CA$subset.cols.y == c(0))			#If NA - set to default
		{
		CA$subset.cols.y <-NULL					#Set to default
		}


						 
	if (!CA$verbose == TRUE)												 
		{
		CA$min.subj <- NULL		 											 
		CA$iterations <- NULL												 
		CA$method.args <- NULL										 		
		CA$verbose <- NULL													 
		CA$iterations.gap <- NULL	
		CA$sim.score.parameters <- NULL	
		CA$subset.cols.x <- NULL
		CA$subset.cols.y <- NULL	
		CA$errthresh <- NULL
		}


	return(CA)
}



method.calculation <-
function(b,nsubj,data,CA){
#*************************************************************************************
#* 	Function to calculate cor  and check that total in cols is not zero              *
#*************************************************************************************
	check.col.sums <- colSums(b)==0
	if (length(check.col.sums[check.col.sums==TRUE]))  		#If there is a column that is all zeros - try to reboot data 5 times
		{
		cnt.tries = 0							#Initialize counter of  tries that we will try to reboot
		while(cnt.tries < CA$retries.max.iterations  && length(check.col.sums[check.col.sums==TRUE]))  #Check if we succeded rebooting the data so no cols are zero
			{
			cnt.tries <- cnt.tries+1			#Increase the counter
			b1 = try.reboot.again (nsubj,data)	#Try to reboot again
			check.col.sums <- colSums(b1)==0	#Are there columns with all zeros in the new rebooted matrix?
			if (!length(check.col.sums[check.col.sums==TRUE])) #If there is no column that is all zeros - post the result to continue processing
				{b <- b1}						#Post the result
			}
		if (cnt.tries > CA$retries.max.iterations  && length(check.col.sums[check.col.sums==TRUE]))	#If reboot did not work - Stop the run
			{
			ErrMsg = paste('Tried to reboot the data', CA$retries.max.iterations, 'times but always get a column that is all zeros')  #Error 
			warning(ErrMsg)

			}
		}
	
	measure.function.parm.list <- append(list(x=b), CA$method.args)	#Build the measure function parameter list
	boot.cor <-do.call(CA$method,measure.function.parm.list)	#Invoke the measuring function		
	return(boot.cor)				
	}

try.reboot.again <-
function(nsubj,data){
#*************************************************************************************
#* 	Function to reboot the data matrix agin so that we don't get cols with all 0     *
#*************************************************************************************
		possible.rows = diag(rep(1,nsubj))
		boot.rowids        = sample(seq(1,nsubj),nsubj,replace=TRUE)
		# Generate the bootstrap matrices by getting the appropriate rows from the possible bootstrap rows
		boot.matrix        = possible.rows[boot.rowids,]
		b1<-resample (data,boot.matrix) 		#b1 is the rebooted matrix that we return
		return(b1)
}
