ccrepe_norm <-
function(data){
#*************************************************************************************
#* A function to normalize a data matrix across the rows                             *
#*************************************************************************************
	data.normalized <- data/apply(data,1,sum)		#Normalize
	data.normalized[is.na(data.normalized )] <- 0 #Replace NA's with Zeros
  	return(data.normalized)
}
