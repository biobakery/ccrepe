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
