merge_two_matrices <-
function(mat1,mat2)
#*************************************************************************************
# A function to merge two matrices, which works when they're different dimensions    *
# merge the two matrices by the rows, merging only rows present in both              *
# (ensures no missing samples)                                                       *
# exclude the first column, which is the row names                                   *
#*************************************************************************************
{
	mat <- as.matrix(merge(mat1,mat2,by="row.names")[,-1])
	class(mat) <- "numeric"    # The default output of merge is strings; convert back to numbers
	return(mat)
}
