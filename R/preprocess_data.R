preprocess_data <-
function(X,SelectedSubset,CA)
#**********************************************************************
#	Preprocess input data 				                              *
#**********************************************************************
{		
		MyDataFrame<-na.omit(X	)									#Post the input data into a working data frame - take only subset requestd by user
		MyDataMatrix<- data.matrix(MyDataFrame, rownames.force = NA)	#Convert into a matrix
		MyDataMatrix1  = MyDataMatrix[,SelectedSubset]					#Select only the columns the User requested (If he did not: subset1=all columns)
		mydata <- MyDataMatrix1[rowSums(MyDataMatrix1 != 0) != 0, ] 	#Remove rows that are all zero to prevent NaNs
		if(nrow(mydata) < CA$min.rows) 						#If not enough data, issue messages in files and stop the run 
			{
			ErrMsg = paste('Not enough data - found ',nrow(mydata),' rows of data - Less rows than  ',CA$min.rows, ' min.rows - Run Stopped')  #Error 
			stop(ErrMsg)
			}
		ProcessedX = mydata/apply(mydata,1,sum)
	return(ProcessedX)
}
