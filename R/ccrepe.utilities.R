calculate_q_values <- function(CA)
#**********************************************************************
#	Calculate Q values    				                              *
#**********************************************************************
{
	q.values.calc.matrix <- CA$p.values   #Set the default p.values  for one and two dataset cases
	if  (CA$OneDataset == TRUE)	
		{q.values.calc.matrix[lower.tri(q.values.calc.matrix)] <- NA}    #If One dataset,  use only the upper part of the matrix

	QValuesArranged <- calculate_q_values_vector(q.values.calc.matrix,CA)$q.values.arranged

	CA$q.values<-q.values.calc.matrix
	CA$q.values[which(!is.na(q.values.calc.matrix))] = QValuesArranged

	return(CA)
}

calculate_q_values_vector <- function(p.values.vector,CA){

	p.values <- p.values.vector
	non.missing.p.values <- p.values[which(!is.na(p.values))]
	m = length(non.missing.p.values)                                #m is the total number of p-values
	ln_m = log(m)									#log m
	ln_m_Plus_Gamma = ln_m + CA$Gamma					
	SortedVector = sort(non.missing.p.values,index.return = TRUE)	#Sort the entire PValues matrix into a vector
	KVector = seq(1,m)						#A vector with the total number of entries in the PValues matrix
	QValues = SortedVector$x*m*ln_m_Plus_Gamma/KVector		#Calculate a vector containing the Q values
	QValuesArranged = rep(-1,m)
	QValuesArranged[SortedVector$ix] = QValues
	q.values <- p.values
	q.values[which(!is.na(p.values))] = QValuesArranged

	return(list(q.values.vec = q.values, q.values.arranged = QValuesArranged))
}


ccrepe_norm <-
function(data){
#*************************************************************************************
#* A function to normalize a data matrix across the rows                             *
#*************************************************************************************
	data.normalized <- data/rowSums(data)		#Normalize
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
	
	CA$p.values <-matrix(data=NA,nrow=n,ncol=n)	#Build the empty PValues matrix
	CA$z.stat <-matrix(data=NA,nrow=n,ncol=n)	#Build the empty z.stat matrix
	
	CA$cor <-matrix(data=NA,nrow=n,ncol=n)	#Build the empty correlation matrix
	
	nsubj = nrow(data)				# Number of subjects

        # The subset of columns for which to calculate the similarity scores
        col.subset <- 1:n
        if( length(CA$subset.cols.1) > 0 && CA$compare.within.x )               #If the User entered a subset of columns
                    {
                    col.subset <- CA$subset.cols.1
                    } else if(length(CA$subset.cols.2) > 0 && !CA$compare.within.x)
                    {
                        col.subset <- unique(c(col.subset,CA$subset.cols.2))
                    }

	# The matrix of possible bootstrap rows (which when multiplied by data give a specific row); of the form with all 0s except for one 1 in each row
	possible.rows = diag(rep(1,nsubj)) 

	#Generating the bootstrap and permutation matrices
	bootstrap.matrices = list()    # The list of matrices of possible bootstrap rows (each matrix multiplies by the data to give a resampled dataset)
	permutation.matrices = list()  # The list of permutation matrices; each matrix will be have columns which are permutations of the row indices
	
	

	if (length(colnames(CA$data1)) == 0)		#If no colum names - force them to be the column number
		{
		colnames(CA$data1)<-1:ncol(CA$data1) 
		}
	
	#**************************************************************
	#*  Pre allocate                                              *
	#**************************************************************
	permutation.matrices  = vector("list", N.rand)		#Pre allocate the permutation matrices
	bootstrap.matrices = vector("list", N.rand)		#Pre allocate

	
	for(i in seq_len(N.rand)){
		if (i %% CA$iterations.gap == 0)   #If output is verbose and the number of iterations is multiple of iterations gap - print status
			{
			print.msg = paste('Sampled ',i,' permutation and bootstrap datasets')
			log.processing.progress(CA,print.msg)  #Log progress
			}
	
   
		# Get the rows of the possible.rows matrix; these correspond to the rows which will be included in the resampled dataset

		boot.rowids = sample(seq(1,nsubj),nsubj,replace=TRUE)

		# Generate the bootstrap matrix by getting the appropriate rows from the possible bootstrap rows
		boot.matrix = possible.rows[boot.rowids,] 

		# Add the bootstrap matrix to the list                   
		bootstrap.matrices [[ i ]] = boot.matrix   #Add bootstrap matrix 

		perm.matrix = replicate(n,sample(seq(1,nsubj),nsubj,replace=FALSE))  # The matrix has each column be a permutation of the row indices
 
		permutation.matrices [[i]] = perm.matrix		#Populate the permutation matrix
	}

	# The bootstrapped data; resample the data using each bootstrap matrix

	
 	log.processing.progress(CA,"Building boot data")  #Log progress
	boot.data = lapply(bootstrap.matrices,resample,data=data)
	

	###  Using method.calculation  for the one dataset also
	
	log.processing.progress(CA,paste("Getting similarity score for bootstrap data: col.subset =",toString(col.subset)))  #Log progress
	boot.cor  = lapply(boot.data,method.calculation,nsubj,data,CA,col.subset )		 ###Function to check is all cols are zeros and apply cor

	log.processing.progress(CA,"Generating permutation data")  #Log progress
	# Generating the permutation data; permute the data using each permutation matrix
	permutation.data = lapply(permutation.matrices,permute,data=data)
	# The correlation matrices for the bootstrapped data; calculate the correlation for each resampled dataset

	# The normalized permutation data; the permutation data needs to be renormalized, but not the bootstrap data
	permutation.norm.raw = lapply(permutation.data,ccrepe_norm)
        permutation.norm     = lapply(permutation.data,function(mat) mat[,col.subset])

	# The correlation matrices of the permuted data; calculate the correlation for each permuted dataset
	
	CA$verbose.requested = FALSE			#If the User requested verbose output - turn it off temporarily
	if (CA$verbose == TRUE)
		{
			CA$verbose.requested <- TRUE		#Turn of the verbose flag save variable
			CA$verbose <- FALSE					#But pass non verbose to the CA$method (Could be nc.score or anything else)	
		}
	
	log.processing.progress(CA,paste("Calculating permutation similarity score: col.subset =",toString(col.subset)))  #Log progress
	permutation.cor <- do.call(lapply,c(list(permutation.norm,CA$method), CA$method.args))  #Invoke the measuring function
	 
	# Now, actually calculating the correlation p-values within the dataset
	log.processing.progress(CA,"Calculating p-values")  #Log progress

	
 
	loop.range1 <- col.subset						#Establish looping range default
	loop.range2 <- col.subset						#Establish looping range default
	max.comparisons <- choose(length(col.subset),2)					#The default number of comparisons
	
	if( length(CA$subset.cols.1) > 0 )		#If the User entered a subset of columns
	    	{
		loop.range1 <- which(col.subset %in% CA$subset.cols.1)		#Use that subset
		

		if( CA$compare.within.x )               #If comparing only within subset.cols.x
		    	{
			loop.range2 <- which(col.subset %in% CA$subset.cols.1)			#Use subset.cols.x for the inner loop as well

			} else if( length(CA$subset.cols.2)>0 ){	   #If comparing between x and y and the user input subset.cols.y

			loop.range2 <- which(col.subset %in% CA$subset.cols.2) #Use subset.cols.y for the inner loop

			}
		}


	if( !is.na(CA$concurrent.output) || CA$make.output.table )
	    {
	    output.table = data.frame(feature1=rep(NA,max.comparisons), 
	    		   	      feature2=rep(NA,max.comparisons), 
				      sim.score=rep(NA,max.comparisons), 
				      z.stat=rep(NA,max.comparisons),
				      p.value=rep(NA,max.comparisons), 
				      q.value=rep(NA,max.comparisons))

	    } else {
	    output.table = NULL
	    }

	if ( !is.na(CA$concurrent.output) )						#If user requested to print the output
	   {
	   cat(colnames(output.table),sep="\t",file=CA$concurrentFile,append=TRUE)
	   cat("\n",file=CA$concurrentFile,append=TRUE)
	   }	

	internal.loop.counter = 0		#Initialize the loop counter
	outer.loop.indices.completed = c()	#Initialize the list to keep track of already completed outer indices

	for( index1 in seq_len(length(loop.range1)) )
	{
		i = loop.range1[index1]
		outer.loop.indices.completed = c(outer.loop.indices.completed,i)	# Keep track of which outer loop indices have been accounted for
		inner.loop.range = setdiff(loop.range2,outer.loop.indices.completed)	# Use as the inner loop only those inner loop indices which haven't been in the outer loop
		{
			for(index2 in seq_len(length(inner.loop.range)))
			{	
				k = inner.loop.range[index2]
				# Get a vector of the (i,k)th element of each correlation matrix in the list of bootstrapped data; this is the bootstrap distribution
			
				internal.loop.counter = internal.loop.counter + 1  #Increment the loop counter
				if (internal.loop.counter %% CA$iterations.gap == 0)   #If output is verbose and the number of iterations is multiple of iterations gap - print status
				{
					print.msg = paste('Completed ', internal.loop.counter, ' comparisons')
					log.processing.progress(CA,print.msg)  #Log progress
				}
			
				 
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

					measure.parameter.list <- append(list(x=data[,col.subset[i]],y=data[,col.subset[k]]), CA$sim.score.parameters)  #build the method do.call parameter list
					cor <- do.call(CA$method,measure.parameter.list)	#Invoke the measuring function

					####################################################
					#  New p.value calculation                         #
					####################################################
					z.stat <- (mean(bootstrap.dist) - mean(permutation.dist))/sqrt(0.5*(var(permutation.dist)+var(bootstrap.dist)))
					p.value <- 2*pnorm(-abs(z.stat))					
					}
				CA$z.stat[col.subset[i],col.subset[k]] = z.stat					#Post z.stat in output matrix	
				CA$z.stat[col.subset[k],col.subset[i]] = z.stat					#Post z.stat in output matrix					
				CA$p.values[col.subset[i],col.subset[k]] = p.value				#Post it in the p-values matrix  
				CA$p.values[col.subset[k],col.subset[i]] = p.value				#Post it in the p-values matrix  
				CA$cor[col.subset[i],col.subset[k]] = cor						#Post it in the cor matrix  
				CA$cor[col.subset[k],col.subset[i]] = cor						#Post it in the cor matrix  				

				
				if( !is.na(CA$concurrent.output) || CA$make.output.table )
				    	   {
					   concurrent.vector <- c(colnames(data)[i],colnames(data)[k],cor,z.stat,p.value,NA)
					   output.table[internal.loop.counter,] = concurrent.vector
					   }

				if ( !is.na(CA$concurrent.output) )						#If user requested to print the output
					   {
				   	   cat(concurrent.vector,sep="\t",file=CA$concurrentFile,append=TRUE)
					   cat("\n",file=CA$concurrentFile,append=TRUE)
					   }	
			}
		}
	}
	

	
	CA <- calculate_q_values(CA)						#Calculate the QValues
	if( CA$make.output.table )
	    {
	    output.table$sim.score <- as.numeric(output.table$sim.score)
	    output.table$z.stat    <- as.numeric(output.table$z.stat)
	    output.table$p.value   <- as.numeric(output.table$p.value)
	    table.q.values <- calculate_q_values_vector(output.table$p.value,CA)$q.values.vec
	    output.table$q.value <- table.q.values
	    }
 
	CA$q.values[lower.tri(CA$q.values)] = t(CA$q.values)[lower.tri(t(CA$q.values))]  #Making the q.values matrix symmetrical for the one dataset case
	
#	for (indx in 1:nrow(data.cor))						#post the q-values
#		{
#			i = data.cor[indx,1]
#			k = data.cor[indx,2]
#			data.cor[indx,5] = CA$q.values[i,k]
#		}
	if( !is.na(CA$concurrent.output) || CA$make.output.table )
	    {
	    CA$output.table <- output.table[1:internal.loop.counter,]		#Post it in the common Area
	    } else 
	    {
	    CA$output.table <- NULL
	    }

	

		
	#********************************************************************
	#*  Final Edits before exiting                                      *
	#********************************************************************
	diag(CA$q.values) <- NA											#Set diagonal of q.values to NA 
	colnames(CA$p.values)<-colnames(CA$data1)						#Set the names of the columns in the p.values matrix
	rownames(CA$p.values)<-colnames(CA$data1)						#Set the names of the rows in the p.values matrix
	colnames(CA$q.values)<-colnames(CA$data1)						#Set the names of the columns in the q.values matrix
	rownames(CA$q.values)<-colnames(CA$data1)						#Set the names of the rows in the q.values matrix
	colnames(CA$cor)<-colnames(CA$data1)							#Set the names of the columns in the q.values matrix
	rownames(CA$cor)<-colnames(CA$data1)							#Set the names of the rows in the q.values matrix
	colnames(CA$z.stat) <- colnames(CA$data1)						#Set the names of the cols in the z.stat matrix
	rownames(CA$z.stat) <- colnames(CA$data1)						#Set the names of the rows in the z.stat matrix

	CA$sim.score <- CA$cor											#Rename cor to sim.score
	diag(CA$p.values) <- NA											#Set diagonal of p.values to NA
        diag(CA$z.stat)   <- NA                                                                                 #Set diagonal of z.stat
	CA$cor <- NULL
	CA <- clean_common_area_after_processing(CA)	#Clean the Common Area before returning to the User
 
	if (length(CA$subset.cols.x) > 1)				#If used a subset - present only the subset
		{
		CA$p.values <- CA$p.values[CA$subset.cols.x,CA$subset.cols.x]   #Display only the subset of cols and rows
		CA$q.values <- CA$q.values[CA$subset.cols.x,CA$subset.cols.x]   #Display only the subset of cols and rows
		CA$sim.score <- CA$sim.score[CA$subset.cols.x,CA$subset.cols.x]   #Display only the subset of cols and rows
                CA$z.stat    <- CA$z.stat[CA$subset.cols.x,CA$subset.cols.x]
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

	log.processing.progress(CA,"Two datasets: Initial merge of the matrices")  #Log progress
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

        # Subset of columns for the bootstrapped dataset
        col.subset <- 1:(n1+n2)
        n1.subset <- n1
	n2.subset <- n2
	if( length(CA$subset.cols.1) > 0 )		#If the User entered a subset of columns
	    	{
		col.subset <- CA$subset.cols.1
                n1.subset  <- length(CA$subset.cols.1)

                if(length(CA$subset.cols.2) > 0)
                       {
                       col.subset <- c(col.subset,n1+CA$subset.cols.2)
		       n2.subset <- length(CA$subset.cols.2)
                       }
                } 

        # Subset of columns for the permutation datasets
        


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
	
	#*****************************************************
	#*  Pre Allocate Matrices                            *
	#*****************************************************
	bootstrap.matrices = vector("list", N.rand)		
	permutation.matrices1 = vector("list", N.rand)	
	permutation.matrices2 = vector("list", N.rand)	

	for(i in seq_len(N.rand)){

		if (i %% CA$iterations.gap == 0)   #If output is verbose and the number of iterations is multiple of iterations gap - print status
			{
			print.msg = paste('Sampled ',i,' permutation and bootstrap datasets')
			log.processing.progress(CA,print.msg)  #Log progress
			}
		
 
		# Get the rows of the two possible.rows matrices; these correspond to the rows which will be included in the resampled datasets
		boot.rowids        = sample(seq(1,nsubj),nsubj,replace=TRUE)

		# Generate the bootstrap matrices by getting the appropriate rows from the possible bootstrap rows
		boot.matrix        = possible.rows[boot.rowids,]

		# Add the bootstrap matrices to the appropriate lists  
			
		bootstrap.matrices[[ i ]] = boot.matrix # Add the bootstrap matrix to the list
		perm.matrix1 = replicate(n1,sample(seq(1,nsubj1),nsubj1,replace=FALSE)) # The matrix has each column be a permutation of the row indices
		permutation.matrices1 [[ i ]] = perm.matrix1      # Add the new matrix to the list
		perm.matrix2 = replicate(n2,sample(seq(1,nsubj2),nsubj2,replace=FALSE)) # The matrix has each column be a permutation of the row indices
		permutation.matrices2 [[ i ]] = perm.matrix2     # Add the new matrix to the list

	}

	# The bootstrapped data; resample the data using each bootstrap matrix
	
	log.processing.progress(CA,"Generating boot data")  #Log progress
	boot.data = lapply(bootstrap.matrices,resample,data=data)
 


	log.processing.progress(CA,"Generating permutation data")  #Log progress
	# Generating the permutation data; permute the data using each permutation matrix
	permutation.data1 = lapply(permutation.matrices1,permute,data=data1)
	permutation.data2 = lapply(permutation.matrices2,permute,data=data2)

	# The normalized permutation data; the permutation data needs to be renormalized, but not the bootstrap data
	permutation.norm1 = lapply(permutation.data1,ccrepe_norm)
	permutation.norm2 = lapply(permutation.data2,ccrepe_norm)


	# The correlation matrices of the permuted and bootstrapped data; see extractCor function for more details
	# mapply is a function that applies over two lists, applying extractCor to the first element of each, then the
	# second element of each, etc.
 
 
	log.processing.progress(CA,paste("Calculating permutation similarity scores: col.subset =",toString(col.subset)))  #Log progress
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
		            endrow=n1.subset,
		            startcol=(n1.subset+1),
		            endcol=length(col.subset), 				 
                	    		MoreArgs = list(col.subset = col.subset,
					my.method=CA$method,
					method.args=CA$method.args,
					outdist=CA$outdist,
					outdistFile=CA$outdistFile),
		            SIMPLIFY=FALSE)
		

	#******************************************************************************** 
	# For each matrix, check if there is a column that is all zeros 				*
	# If such matrix is found,  try to reboot it at most 5 times until such matrix  *
	# is found that does not contain a column with all zeros                        *
	# The limit of 5 tries is now hard coded and needs to be re-evaluated           *
	# If after 5 times there is no success - we stop the run (Need to verify!!!!    *
	#********************************************************************************
	
	log.processing.progress(CA,paste("Calculating boot similarity scores: col.subset  =",toString(col.subset)))  #Log progress
	boot.cor  = lapply(boot.data,method.calculation,nsubj,data,CA,col.subset )		 ###Function to check is all cols are zeros and apply cor
	
	# Now calculating the correlations and p-values between the two datasets
	log.processing.progress(CA,"Calculating p-values")  #Log progress

	
	
	CA$p.values <-matrix(data=NA,nrow=n1,ncol=n2)	#Build the empty p.values matrix
	CA$z.stat  <-matrix(data=NA,nrow=n1,ncol=n2)		#Build the empty z.stat matrix
	CA$cor <-matrix(data=NA,nrow=n1,ncol=n2)	#Build the empty correlation matrix
	
	loop.range1 <- 1:n1						#Establish looping range default
	
	if ( length(CA$subset.cols.1)> 0 )		#If the User entered a subset of columns
		{
		loop.range1 <- which(col.subset %in% CA$subset.cols.1)		#Use the subset of columns
		}

	loop.range2 <- 1:n2						#Establish looping range default
	
	if ( length(CA$subset.cols.2) > 0 )		#If the User entered a subset of columns
		{
		loop.range2 <- which(col.subset %in% (n1 + CA$subset.cols.2)) - n1.subset		#Use the subset of columns
		}


	max.comparisons <- n1*n2
	if( !is.na(CA$concurrent.output) || CA$make.output.table )
	    {
	    output.table = data.frame(feature1=rep(NA,max.comparisons), 
		       		      feature2=rep(NA,max.comparisons), 
				      sim.score=rep(NA,max.comparisons), 
				      z.stat=rep(NA,max.comparisons),
				      p.value=rep(NA,max.comparisons), 
				      q.value=rep(NA,max.comparisons))
	    }

	if ( !is.na(CA$concurrent.output) )						#If user requested to print the output
	   {
	   cat(colnames(output.table),sep="\t",file=CA$concurrentFile,append=TRUE)
	   cat("\n",file=CA$concurrentFile,append=TRUE)
	   }	

	internal.loop.counter = 0   # Initialize the counter
	for(index1 in seq_along(loop.range1))
	{
		i = loop.range1[index1]
		for(index2 in seq_along(loop.range2))
		{
			k = loop.range2[index2]
			# Get a vector of the (i,k)th element of each correlation matrix in the list of bootstrapped data; this is the bootstrap distribution

				internal.loop.counter = internal.loop.counter + 1  #Increment the loop counter
				if (internal.loop.counter %% CA$iterations.gap == 0)   #If output is verbose and the number of iterations is multiple of iter
#ations gap - print status
				{
					print.msg = paste('Completed ', internal.loop.counter, ' comparisons')
					log.processing.progress(CA,print.msg)  #Log progress
				}
				
			bootstrap.dist = unlist(lapply(boot.cor,'[',i,n1.subset+k))
			bootstrap.dist[is.na(bootstrap.dist)] <- 0				#If there is an NA in bootstrap.dist - replace with 0 (Needs 
#review)
			
			# Get a vector the (i,k)th element of each correlation matrix in the list of permuted data; this is the permuted distribution
			permutation.dist = unlist(lapply(permutation.cor,'[',i,k))
			
			if    (!is.na(CA$outdist))						#If user requested to print the distributions
					{
					RC <- print.dist(bootstrap.dist,permutation.dist,CA,i,k)
					}

			n.0_1 = sum(data[,col.subset[i]]==0)				#Number of zeros in column i
			n.0_2 = sum(data[,(col.subset[n1.subset+k]-n1)]==0)				#Number of zeros in column n1+k
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
						measure.parameter.list <- append(list(x=data[,col.subset[i]],y=data[,col.subset[n1.subset+k]]), CA$sim.score.parameters)  #build the method 
#do.call parameter list
						cor <- do.call(CA$method,measure.parameter.list)	#Invoke the measuring function

						####################################################
						#  New p.value calculation                         #
						####################################################
					z.stat <- (mean(bootstrap.dist) - mean(permutation.dist))/sqrt(0.5*(var(permutation.dist)+var(bootstrap.dist)))
						p.value <- 2*pnorm(-abs(z.stat))					
					}

								  
			CA$p.values[col.subset[i],col.subset[n1.subset+k]-n1] = p.value				#Post it in the p.values matrix  
			CA$z.stat[col.subset[i],col.subset[n1.subset+k]-n1] = z.stat					#Post it in the z.stat matrix 
			CA$cor[col.subset[i],col.subset[n1.subset+k]-n1] = cor					#Post it in the cor matrix  

			if( !is.na(CA$concurrent.output) || CA$make.output.table )
			    {
			    concurrent.vector <- c(colnames(data1)[col.subset[i]],colnames(data2)[col.subset[n1.subset+k]-n1],cor,z.stat,p.value,NA)	  
			    output.table[internal.loop.counter,] = concurrent.vector
			    }

			if ( !is.na(CA$concurrent.output) )						#If user requested to print the output
			   {
		   	   cat(concurrent.vector,sep="\t",file=CA$concurrentFile,append=TRUE)
			   cat("\n",file=CA$concurrentFile,append=TRUE)
			   }	

			

		}
	}
	


	CA <- calculate_q_values(CA)						#Calculate the QValues

	if( CA$make.output.table )
	    {
	    output.table$sim.score <- as.numeric(output.table$sim.score)
	    output.table$z.stat    <- as.numeric(output.table$z.stat)
	    output.table$p.value   <- as.numeric(output.table$p.value)
	    table.q.values <- calculate_q_values_vector(output.table$p.value,CA)$q.values.vec
	    output.table$q.value <- table.q.values
	    }

	
#	for (indx in 1:nrow(CA$data.cor))						#post the q-values
#		{
#			i = CA$data.cor[indx,1]
#			k = CA$data.cor[indx,2]
#			CA$data.cor[indx,5] = CA$q.values[i,k]
#		}

	if( !is.na(CA$concurrent.output) || CA$make.output.table )
	    {
	    CA$output.table <- output.table[1:internal.loop.counter,]
	    } else
	    {
	    CA$output.table <- NULL
	    }
	
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
        CA$z.stat    <- CA$z.stat[total.rows.to.display,total.cols.to.display]
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
		CA$subset.cols.1 = NULL					#Set to default
		}
	if  ( length(CA$subset.cols.2) == 0 && (is.null(CA$subset.cols.2)))				#If NULL - set to default
		{
		CA$subset.cols.2 = NULL					#Set to default
		}
	
  
	if    (!is.na(CA$outdist))							#If the user passed a file - open it
		{
		CA$outdistFile = file(CA$outdist,open='at')		#Open outdist file
		}

	if    (!is.na(CA$concurrent.output))							#If the user passed a file - open it
		{
		CA$concurrentFile = file(CA$concurrent.output,open='at')		#Open the concurrent output file
		}

	if ( CA$compare.within.x !=  TRUE  & 	CA$compare.within.x !=  FALSE )	#compare.within.x must be either true or false
		{
		CA$compare.within.x =  TRUE									#True - is the default
		}

	if ( CA$make.output.table !=  TRUE  & 	CA$make.output.table !=  FALSE )	#make.output.table flag must be either true or false
		{
		CA$make.output.table =  FALSE									#False - is the defauls
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
	CA$sim.score.parameters<-list()							#Initialize the parameter list
	for (name in names(CA$method.args)) {					#Add the entries in method.args to the measuring parameter list			
 		CA$sim.score.parameters[[name]]<-CA$method.args[[name]]
		}	
		

	CA$retries.max.iterations =  -round(log2(CA$errthresh)) #This is the maximum number of iterations to try to reboot a matrix if in a col all values are 0
		
		
	return(CA)			 				#Return list of decoded input parameters
}




extractCor <-
function(mat1,mat2,startrow,endrow,startcol,endcol,col.subset,my.method,method.args,outdist,outdistFile,  ...)
#******************************************************************************************
# A function to calculate the correlation of the two matrices by merging them,            *
#     calculating the correlation of the merged matrix, and extracting the appropriate    *
#     submatrix to obtain the correlation of interest                                     *
# startrow, endrow, startcol, endcol give the indeces of the submatrix to extract         *
# The method and the method parameters are passed via the list
#******************************************************************************************
{
  	mat <- merge_two_matrices(mat1,mat2)[,col.subset]	            #Merge the two matrices
	measure.function.parm.list <- append(list(x=mat), method.args)	
	mat_C <-do.call(my.method,measure.function.parm.list)	#Invoke the measuring fnction
        if(length(startrow:endrow)==1){
            sub_mat_C <- matrix(mat_C[startrow:endrow,startcol:endcol],nrow=1)
        } else if(length(startcol:endcol)==1){
            sub_mat_C <- matrix(mat_C[startrow:endrow,startcol:endcol],ncol=1)
        } else {
            sub_mat_C <- mat_C[startrow:endrow, startcol:endcol] # Extract the appropriate submatrix
        }
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
		ProcessedX = mydata/rowSums(mydata)
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

	if   (!is.na(CA$concurrent.output))										#If user requested to print the distributions
		{
		close(CA$concurrentFile)										#Close outdist file	
		CA$concurrentFile <- NULL										#And remove it from the common area
		}

	if( !CA$make.output.table )
	    {
	    CA$output.table <- NULL
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
	CA$concurrent.output <- NULL
	CA$Gamma <- NULL 
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
		CA$compare.within.x <- NULL
		CA$sim.score.parameters <- NULL	
		CA$subset.cols.x <- NULL
		CA$subset.cols.y <- NULL	
		CA$errthresh <- NULL
		CA$make.output.table <- NULL
		}



	return(CA)
}



method.calculation <-
function(b,nsubj,data,CA,col.subset){
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
	
	measure.function.parm.list <- append(list(x=b[,col.subset]), CA$method.args)	#Build the measure function parameter list
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

log.processing.progress <-
function(CA,msg){
#*************************************************************************************
#* 	Function to report the progress of processing                                    *
#*  only if verbose flag is on                                                       *
#*************************************************************************************
		if (CA$verbose == TRUE  ||  (!is.null(CA$verbose.requested)  && CA$verbose.requested == TRUE ))
			message(cat(date(),' ==> ',msg))			#Display date and time and progress
		return (0)
}
