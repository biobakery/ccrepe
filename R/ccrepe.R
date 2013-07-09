ccrepe <-
function(
#*************************************************************************************
#*  	ccrepe                                                                       *
#*************************************************************************************
	x=NA,								#Data Frame  1 - For ccrepe and nc.score
	y=NA,								#Data Frame  2 - For ccrepe and nc.score 
	method = cor,						#Default
	method.args = list(), 				#Arguments for the method
	min.subj = 20,						#Minimum rows in "data" frame to proceed (If not - run stops) - For ccrepe
	iterations = 1000,					#Reboot iterations - For ccrepe
	subset.cols.1 = c(0),				#Subset of cols from cav1 to iterate on (c(0)== ALL) - For ccrepe
	subset.cols.2 = c(0),				#Subset of cols from cav2 to iterate on (c(0)== ALL) - For ccrepe
	errthresh = 0.0001, 				#Threshold error if there is enough data to calculate cor an pval for certain i and k - For ccrepe
	verbose = FALSE,					#Request for verbose output
	iterations.gap = 100,				#If output is verbose - after how many iterations issue a status message?
	distributions = NA					#Output Distribution file - For ccrepe
	)					

{
	#**********************************************************************
	# Invoke ccrepe using the method                                      *
	#**********************************************************************
	

	CA <-list(data1=x,
			data2=y,
			min.subj=min.subj,
			iterations=iterations,
			subset.cols.1=subset.cols.1,
			subset.cols.2=subset.cols.2,
			errthresh=errthresh,
			method=method,
			method.args=method.args,
			verbose=verbose,
			iterations.gap=iterations.gap,
			outdist=distributions)

			
	CA <- ccrepe.calculations(CA)
	return(CA)					

}




