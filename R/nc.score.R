nc.score <-
function(
#*************************************************************************************
#*  	NC_Score - Start, invoke the function, return  Output values                 *
#*************************************************************************************
    frmOne=NA,							#First input 
	frmTwo=NA,						    #Second input 
	adBins = NA,
	min.abundance = 0.0001,
	min.samples = 0.1
		)
	{
	CA <- preprocess_nc_score_input(frmOne,             #Preprocess input and build the common area
			frmTwo,
			adBins,
			min.abundance,
			min.samples)
#*************************************************************************************
#*  	NC_Score                                                                     *
#*************************************************************************************
  data <- CA$data1_trimmed_cat_matrix			#Post the input data to be processed
  CA$Output$NCScoreDetail <- data.frame()		#Define Output List as an empty dataframe				
  CA$Output$NCScore.matrix <-matrix(nrow=nrow(data),ncol=nrow(data))	#Define output NCScore matrix
  mode(data) <- "numeric"
 
  n <-  CA$adBins								#n is the number of bins

  adj <- ((1.5)*n*(n-1)/(n^2-n+1))
  for(i in 1:(nrow(data))) {
    for(j in (i):nrow(data)) {
      ijsum <- 0
      cosum <- 0
      cesum <- 0
      for (m in 1:(ncol(data)-1)) {
        for (n in (m+1):ncol(data)) {
          mx <- (c(data[i,m], data[j,m], data[i,n], data[j,n]))
          if (length(unique(mx)) >= 2) {
            if (((mx[1]<mx[3])&(mx[2]<mx[4])) | ((mx[1]>mx[3])&(mx[2]>mx[4]))) {
              cosum <- cosum + 1 }
            if (((mx[1]>mx[3])&(mx[3]<mx[4])&(mx[4]>mx[2])&(mx[2]<mx[1])) | ((mx[1]<mx[3])&(mx[3]>mx[4])&(mx[4]<mx[2])&(mx[2]>mx[1]))) {
              cesum <- cesum + 1 }  
          }
        }
      }

	  
      cesum_adj <- cesum * adj
      ijsum <- (cosum - cesum_adj)
      CA$Output$NCScore.matrix[i,j] <- ijsum				#Post the result
      CA$Output$NCScore.matrix[j,i] <- ijsum				#Post the result in the symmetrical side
      
	  if (is.null(rownames(data)))							#Check if row name is null
		  {rownames(data)[1:nrow(data)] <- paste("Row", 1:nrow(data), sep="")}
		          	   
		  
      output.result <-data.frame(list(
		i=i,
		j=j,		
		Bug1=rownames(data)[i],
		Bug2=rownames(data)[j],
		IJSum = ijsum))

      CA$Output$NCScoreDetail <- rbind(CA$Output$NCScoreDetail,output.result)	#Build an output row to post the results	
               
    }
  }
 
  NS	<- ncol(data)													#Number of Samples
  NB 	<-   CA$adBins													#Number of bins
  RenormalizationFactor <-  choose(NS, 2) - (NS %% NB) * choose((floor(NS/NB) + 1), 2) - (NB - NS %% NB) * choose(floor(NS/NB), 2)
  if  (RenormalizationFactor == 0) 										#So that we dont get NaNs
		{RenormalizationFactor =  1} 

  CA$Output$NCScore.matrix <-   CA$Output$NCScore.matrix / RenormalizationFactor  #Renormalize the matrix
  if (CA$YEntered == TRUE) {										#Vector Calculation
		return (CA$Output$NCScore.matrix[1,2])						#The response is a Number
		}
  return(CA$Output$NCScore.matrix)
	}
