ccrepe_process_one_dataset <-
function(data,N.rand, CA){
#*************************************************************************************
#* 	Function to apply the new method to a dataset                                    *
#*  As a note, data needs to be a matrix with no missing values                      *
#*************************************************************************************

	time_start = Sys.time()
	n = ncol(data)					# Number of columns, starting at 1; this is also the number of bugs
	
	CA$p.values <-matrix(data=0,nrow=n,ncol=n)	#Build the empty PValues matrix
	
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
	boot.cor = lapply(boot.data,cor,method=CA$method)


	# Generating the permutation data; permute the data using each permutation matrix
	permutation.data = lapply(permutation.matrices,permute,data=data)
	# The correlation matrices for the bootstrapped data; calculate the correlation for each resampled dataset

	# The normalized permutation data; the permutation data needs to be renormalized, but not the bootstrap data
	permutation.norm = lapply(permutation.data,ccrepe_norm)

	# The correlation matrices of the permuted data; calculate the correlation for each permuted dataset

	permutation.cor = lapply(permutation.norm,cor,method=CA$method)

	
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
					cor = cor(data[,i],data[,k], method=CA$method)   # Calculate the  correlation between the bugs
					# The Z-test to get the p-value for this comparison; as in get.renorm.null.pval
					p.value = pnorm(mean(permutation.dist), 
                                        mean=mean(bootstrap.dist), 
                                        sd=sqrt((var(bootstrap.dist) + var(permutation.dist))*0.5))
					}
				CA$p.values[i,k] = p.value				#Post it in the p-values matrix  
				CA$p.values[k,i] = p.value				#Post it in the p-values matrix  
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
		
	diag(CA$p.values) <- NA											#Set diagonal of p.values to NA 

	time_end = Sys.time()
	return(CA)														# Return the output matrix
}
