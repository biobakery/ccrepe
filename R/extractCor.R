extractCor <-
function(mat1,mat2,startrow,endrow,startcol,endcol,method)
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
