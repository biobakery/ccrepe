ccrepe <-
function(
#*************************************************************************************
#*  	ccrepe                                                                       *
#*************************************************************************************
	x=NA,								#Data Frame  1 - For ccrepe and nc.score
	y=NA,								#Data Frame  2 - For ccrepe and nc.score 
	method = cor,						#Default
	method.args = list(), 	 			#Could be, for example :   list(	method = "spearman")
	min.subj = 20,						#Minimum rows in "data" frame to proceed (If not - run stops) - For ccrepe
	iterations = 1000,					#Reboot iterations - For ccrepe
	subset.cols.1 = c(0),				#Subset of cols from cav1 to iterate on (c(0)== ALL) - For ccrepe
	subset.cols.2 = c(0),				#Subset of cols from cav2 to iterate on (c(0)== ALL) - For ccrepe
	errthresh = 0.0001, 				#Threshold error if there is enough data to calculate cor an pval for certain i and k - For ccrepe
	verbose = FALSE,					#Request for verbose output?
	iterations.gap = 100,				#If output is verbose - after how many iterations issue a status message?
	distributions = NA					#Output Distribution file - For ccrepe
	)					

{
	#**********************************************************************
	# Invoke ccrepe using the method                                      *
	#**********************************************************************
			CA <- ccrepe.calculations(
				x, 							#Data Frame  1 - For ccrepe and nc.score
				y,							#Data Frame  2 - For ccrepe and nc.score 
				method,						#cor 
				method.args$method,			#correlation method
				min.subj,					#Minimum rows in "data" frame to proceed (If not - run stops) - For ccrepe
				iterations,					#Reboot iterations - For ccrepe
				subset.cols.1,				#Subset of cols from cav1 to iterate on (c(0)== ALL) - For ccrepe
				subset.cols.2,				#Subset of cols from cav2 to iterate on (c(0)== ALL) - For ccrepe
				errthresh, 					#Threshold error if there is enough data to calculate cor an pval for certain i and k - For ccrepe
				verbose,					#Request for verbose output?
				iterations.gap,				#If output is verbose - after how many iterations issue a status message?
				distributions  				#Output Distribution file - For ccrepe
				)
				return(CA)					
		
		

}




ccrepe.calculations  <-
#*************************************************************************************
#*  	ccrepe calculations                                                          *
#*************************************************************************************
function(
 				x, 							#Data Frame  1 - For ccrepe and nc.score
				y,							#Data Frame  2 - For ccrepe and nc.score 
				method,						#Parameters,
				method.args.method,			#correlation method
				min.subj,					#Minimum rows in "data" frame to proceed (If not - run stops) - For ccrepe
				iterations,					#Reboot iterations - For ccrepe
				subset.cols.1,				#Subset of cols from cav1 to iterate on (c(0)== ALL) - For ccrepe
				subset.cols.2,				#Subset of cols from cav2 to iterate on (c(0)== ALL) - For ccrepe
				errthresh, 					#Threshold error if there is enough data to calculate cor an pval for certain i and k - For ccrepe
				verbose,					#Request for verbose output?
				iterations.gap,				#If output is verbose - after how many iterations issue a status message?
				distributions 				#Output Distribution file - For ccrepe
		)
{
	Output <-list()
	CA <-list(	data1=x,
			data2=y,
			min.subj=min.subj,
			iterations=iterations,
			errthresh=errthresh,
			method=method,
			method.args.method=method.args.method,
			subset.cols.1=subset.cols.1,
			subset.cols.2=subset.cols.2,
			verbose=verbose,
			iterations.gap=iterations.gap,
			outdist=distributions
			)

	
	CA <- DecodeInputCAArguments(CA)							#Decode Input Parameters

	if (CA$OneDataset == TRUE)
		{
		mydata.norm = preprocess_data(CA$data1,CA$subset.cols.1,CA)						#Preprocess the data 
		CA  = ccrepe_process_one_dataset(mydata.norm,CA$iterations, CA)  				#Invoke the new function 
		CA$data1 <- NULL										#Remove it from the output
		CA$data2 <- NULL										#Remove it from the output
		CA$subset.cols.2<-NULL									#Remove it from output
		CA$outdist <- NULL										#Remove it from output
		CA$OneDataset <- NULL									#Remove it from output
		CA$method <- NULL									#Remove it from output
		}
	else
		{
		CA$data1.norm = preprocess_data(CA$data1,CA$subset.cols.1,CA)					#Preprocess data1 
		CA$data2.norm = preprocess_data(CA$data2,CA$subset.cols.2,CA)					#Preprocess data2
		CA = ccrepe_process_two_datasets  (CA$data1.norm ,CA$data2.norm ,CA$iterations, CA)
		}
	return ( CA )
}






