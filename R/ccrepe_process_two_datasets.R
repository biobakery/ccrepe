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
	boot.cor = lapply(boot.data,cor,use="complete.obs",method="spearman")


	# Now calculating the correlations and p-values between the two datasets
    n.c = 0	# Counter for the number of comparisons (to enter in the output matrix)
	for(i in 1:n1){
		for(k in 1:n2){
			# Get a vector of the (i,k)th element of each correlation matrix in the list of bootstrapped data; this is the bootstrap distribution
			bootstrap.dist = unlist(lapply(boot.cor,'[',i,n1+k))

			# Get a vector the (i,k)th element of each correlation matrix in the list of permuted data; this is the permuted distribution
			permutation.dist = unlist(lapply(permutation.cor,'[',i,k))

			n.c = n.c + 1
			cor.meas[n.c] = cor(data[,i],data[,n1+k],use="complete.obs",method="spearman")   # Calculate the Spearman correlation between the bugs

			# The Z-test to get the p-value for this comparison; as in get.renorm.null.pval
			p.value = pnorm(mean(permutation.dist), 
                                  mean=mean(bootstrap.dist), 
                                  sd=sqrt((var(bootstrap.dist) + var(permutation.dist))*0.5))

			p.values[n.c] = p.value
			bug1[n.c] = colnames(data1)[i]
			bug2[n.c] = colnames(data2)[k]

		}
	}
	data.cor <- data.frame(bug1,bug2,cor.meas,p.values)
	return(data.cor)			# Return the output matrix
}
