qc_filter <-
function(data,CA) {
#********************************************************************************************
#*  	quality control: select species with >= 1E-4 relative abundance in >= 10% of samples*
#********************************************************************************************
	tmp <- {}
	names <- {}

	for (i in 1:ncol(data)) 
	{
		if (length(which(data[,i] >= CA$min.abundance)) >= CA$min.samples*nrow(data))
		{
			tmp <- cbind(tmp, data[,i])
			names <- c(names, colnames(data)[i])
		}
	}
  colnames(tmp) <- names					#Post the column names
  rownames(tmp) <- rownames(data)			#Post the data row names
  return(tmp)
}
