ccrepe_process_two_datasets <-
function(data1,data2,N.rand, CA){
#*************************************************************************************
#* 	Function to apply the new method to two datasets                                 *
#*************************************************************************************
	# Expanding the new method to run on two datasets; this was more complicated than expected.

	time_start = Sys.time()

	# Get number of bugs, subjects for each dataset
	n1 = ncol(data1)
	nsubj1 = nrow(data1)
	n2 = ncol(data2)
	nsubj2 = nrow(data2)

	# Merge the data (so we can calculate the correlations appropriately for each pair of bugs)
	data = merge_two_matrices(data1,data2)

	# Generating the output data matrix (I chose this form for simplicity; we do want to keep the output object George has made)
    data.cor = matrix(ncol=4,nrow=n1*n2)
	colnames(data.cor) = c("bug1", "bug2", "cor", "p.value")

	# The matrix of possible bootstrap rows (which when multiplied by data give a specific row); of the form with all 0s except for one 1 in each row
	possible.rows1 = diag(rep(1,nsubj1)) 
	possible.rows2 = diag(rep(1,nsubj2))

	#Generating the bootstrap and permutation matrices
	bootstrap.matrices1 = list()    # The list of matrices of possible bootstrap rows (each matrix multiplies by the data to give a resampled dataset)
	permutation.matrices1 = list()  # The list of permutation matrices; each matrix will be have columns which are permutations of the row indices
	bootstrap.matrices2 = list()    # The list of matrices of possible bootstrap rows (each matrix multiplies by the data to give a resampled dataset)
	permutation.matrices2 = list()  # The list of permutation matrices; each matrix will be have columns which are permutations of the row indices

	for(i in 1:N.rand){

		# Get the rows of the two possible.rows matrices; these correspond to the rows which will be included in the resampled datasets
		boot.rowids1        = sample(seq(1,nsubj1),nsubj1,replace=TRUE)
		boot.rowids2        = sample(seq(1,nsubj2),nsubj2,replace=TRUE)

		# Generate the bootstrap matrices by getting the appropriate rows from the possible bootstrap rows
		boot.matrix1        = possible.rows1[boot.rowids1,]
		boot.matrix2        = possible.rows2[boot.rowids2,] 

		# Add the bootstrap matrices to the appropriate lists  
		bootstrap.matrices1 = lappend(bootstrap.matrices1,boot.matrix1) # Add the bootstrap matrix to the list
		bootstrap.matrices2 = lappend(bootstrap.matrices2,boot.matrix2) # Add the bootstrap matrix to the list

		perm.matrix1 = replicate(n1,sample(seq(1,nsubj1),nsubj1,replace=FALSE)) # The matrix has each column be a permutation of the row indices
		permutation.matrices1 = lappend(permutation.matrices1,perm.matrix1)     # Add the new matrix to the list
		perm.matrix2 = replicate(n2,sample(seq(1,nsubj2),nsubj2,replace=FALSE)) # The matrix has each column be a permutation of the row indices
		permutation.matrices2 = lappend(permutation.matrices2,perm.matrix2)     # Add the new matrix to the list
	}

	
	# The bootstrapped data; resample the data using each bootstrap matrix
	boot.data1 = lapply(bootstrap.matrices1,resample,data=data1)
	boot.data2 = lapply(bootstrap.matrices2,resample,data=data2)

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
		                   method="spearman",
		                   SIMPLIFY=FALSE)
	boot.cor = mapply(extractCor,
		                   mat1=boot.data1,
		                   mat2=boot.data2,
		                   startrow=1,
		                   endrow=n1,
		                   startcol=(n1+1),
		                   endcol=(n1+n2),
		                   method="spearman",
		                   SIMPLIFY=FALSE)

	# Now calculating the correlations and p-values between the two datasets
    n.c = 0	# Counter for the number of comparisons (to enter in the output matrix)
	for(i in 1:n1){
		for(k in 1:n2){
			# Get a vector of the (i,k)th element of each correlation matrix in the list of bootstrapped data; this is the bootstrap distribution
			bootstrap.dist = unlist(lapply(boot.cor,'[',i,k))

			# Get a vector the (i,k)th element of each correlation matrix in the list of permuted data; this is the permuted distribution
			permutation.dist = unlist(lapply(permutation.cor,'[',i,k))

			cor = cor(data[,i],data[,n1+k],use="complete.obs",method="spearman")   # Calculate the Spearman correlation between the bugs

			# The Z-test to get the p-value for this comparison; as in get.renorm.null.pval
			p.value = pnorm(mean(permutation.dist), 
                                  mean=mean(bootstrap.dist), 
                                  sd=sqrt((var(bootstrap.dist) + var(permutation.dist))*0.5))

			n.c = n.c + 1
			data.cor[n.c,] = c(colnames(data1)[i],colnames(data2)[k],cor,p.value)
		}
	}
	CA$data.cor = data.cor			# Post the result


	time_end = Sys.time()
	return(CA)						# Return the output  
}
