ccrepe <-
function(
#*************************************************************************************
#*  	ccrepe                                                                       *
#*************************************************************************************
	x=NA,								#Data Frame  1 - For ccrepe and nc.score
	y=NA,								#Data Frame  2 - For ccrepe and nc.score 
	sim.score = cor,						#Default
	sim.score.args = list(), 				#Arguments for the method
	min.subj = 20,						#Minimum rows in "data" frame to proceed (If not - run stops) - For ccrepe
	iterations = 1000,					#Reboot iterations - For ccrepe
	subset.cols.x = NULL,				#Subset of cols from cav1 to iterate on (NULL = ALL) - For ccrepe
	subset.cols.y = NULL,				#Subset of cols from cav2 to iterate on (NULL = ALL) - For ccrepe
	errthresh = 0.0001, 				#Threshold error if there is enough data to calculate cor an pval for certain i and k - For first dataset
	verbose = FALSE,					#Request for verbose output
	iterations.gap = 100,				#If output is verbose - after how many iterations issue a status message?
	distributions = NA,					#Output Distribution file - For ccrepe
	compare.within.x=TRUE,				#Boolean of whether to do comparisons within subset.cols.x or between subset.cols.x and subset.cols.y
	concurrent.output = NA,				#If specified, file for concurrent output of the table
	make.output.table = FALSE			#Whether to have the output in table format
	)					

{
	#**********************************************************************
	# Invoke ccrepe using the method                                      *
	#**********************************************************************
	

	CA <-list(data1=x,
			data2=y,
			min.subj=min.subj,
			iterations=iterations,
			subset.cols.1=subset.cols.x,
			subset.cols.2=subset.cols.y,
			compare.within.x=compare.within.x,
			errthresh=errthresh,
			method=sim.score,
			method.args=sim.score.args,
			verbose=verbose,
			iterations.gap=iterations.gap,
			outdist=distributions,
			concurrent.output=concurrent.output,
			make.output.table=make.output.table)

			
	CA <- ccrepe.calculations(CA)
	return(CA)					

}




