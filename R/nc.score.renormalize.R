#!/usr/bin/env Rscript
nc.score.renormalize <-
function(
#*************************************************************************************
#*  	nc.score.renormalize                                                         *
#*************************************************************************************
   	x,						#First discretized entity 
	y,						#Second discretized entity
	nc.score.result)					  
{
	NS	<- nrow(as.matrix(x))										#Number of Samples
	NB 	<-   length(unique(c(x,y)))									#Number of bins
	renormalization.factor <-  choose(NS, 2) - (NS %% NB) * choose((floor(NS/NB) + 1), 2) - (NB - NS %% NB) * choose(floor(NS/NB), 2)
	if  (renormalization.factor == 0) 								#So that we dont get NaNs
		{renormalization.factor =  1} 		
	nc.score.result <- nc.score.result  / renormalization.factor	#Normalize
	return(nc.score.result)
}