qc_filter <-
function(data,CA) {
#********************************************************************************************
#*  	quality control: select species with >= 1E-4 relative abundance in >= 10% of samples*
#********************************************************************************************
  tmp <- {}
  names <- {}

  for (i in 1:nrow(data)) {
	if (length(which(data[i,] >= CA$min.abundance)) >= CA$min.samples*ncol(data)) {
      tmp <- rbind(tmp, data[i,])
      names <- c(names, rownames(data)[i])
    }
  }
  rownames(tmp) <- names
  return(tmp)
}
