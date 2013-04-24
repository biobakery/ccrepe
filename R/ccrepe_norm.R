ccrepe_norm <-
function(data){
#*************************************************************************************
#* A function to normalize a data matrix across the rows                             *
#*************************************************************************************
	return(data/apply(data,1,sum))
}
