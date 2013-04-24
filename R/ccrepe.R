ccrepe <-
function(
#*************************************************************************************
#*  	ccrepe main function                                                         *
#*************************************************************************************
	frmOne=NA,							#Data Frame  1 - For ccrepe and nc.score
	frmTwo=NA,							#Data Frame  2 - For ccrepe and nc.score 
	method = 'pearson',					#Type of correlation - For ccrepe
	min.rows = 20,						#Minimum rows in "data" frame to proceed (If not - run stops) - For ccrepe
	Iterations = 1000,					#Reboot Iterations - For ccrepe
	subset1 = c(0),						#Subset of cols from cav1 to iterate on (c(0)== ALL) - For ccrepe
	subset2 = c(0),						#Subset of cols from cav2 to iterate on (c(0)== ALL) - For ccrepe
	errthresh = 0.0001, 				#Threshold error if there is enough data to calculate cor an pval for certain i and k - For ccrepe
	distributions = NA					#Output Distribution file - For ccrepe
	)

	{

	#**********************************************************************
	#Build a list called "CA" containing input parameters                 *
	#**********************************************************************
	Output <-list()
	CA <-list(	data1=frmOne,
			data2=frmTwo,
			min.rows=min.rows,
			Iterations=Iterations,
			errthresh=errthresh,
			method=method,
			subset1=subset1,
			subset2=subset2,
			Output=Output,
			outdist=distributions
			)

	
	CA <- DecodeInputCAArguments(CA)							#Decode Input Parameters

	if (CA$OneDataset == TRUE)
		{
		mydata.norm = preprocess_data(CA$data1,CA$subset1,CA)							#Preprocess the data 
		CA  = ccrepe_process_one_dataset(mydata.norm,CA$Iterations, CA)  				#Invoke the new function 
		}
	else
		{
		CA$data1.norm = preprocess_data(CA$data1,CA$subset1,CA)					#Preprocess data1 
		CA$data2.norm = preprocess_data(CA$data2,CA$subset2,CA)					#Preprocess data2
		CA = ccrepe_process_two_datasets  (CA$data1.norm ,CA$data2.norm ,CA$Iterations, CA)
		}

	return ( CA$data.cor )
	}
