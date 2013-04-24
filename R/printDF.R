printDF <-
function(Input1, CA, DistributionType)	{
#*************************************************************************************
#*  Function to print a data frame                                                   *
#*************************************************************************************
	OutputLine = c('Distribution:',DistributionType  )
	cat(OutputLine,file=CA$outdistFile,append=T)		#Output header to outdist file
	cat('\n',file=CA$outdistFile,append=T)				#Line feed
	df <- data.frame(Input1)							#convert input to data frame
	pm <- data.matrix(df, rownames.force = NA)			#convert df to matrix
	write.table(pm,CA$outdistFile,quote=F,row.names=F,col.names=F) 	#Write to file
	cat('\n',file=CA$outdistFile,append=T)		   		#Line feed
	return (0)
	}
