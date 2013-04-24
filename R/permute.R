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
