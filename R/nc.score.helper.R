nc.score.helper <-
function(
#*************************************************************************************
#*  	nc.score.helper                                                              *
#*************************************************************************************
   	x,						#First discretized  input 
	y)						#Second discretized input 
{
	ijsum <- 0				#Reset ijsum
	cosum <- 0				#Reset cosum
	cesum <- 0				#Reset cesum
	n <- length(unique(c(x)))	#<------ Need to QA and verify this!!! GW
	adj <- ((1.5)*n*(n-1)/(n^2-n+1))
	for (i in 1:(length(x)-1)) 		#Loop on the entries of x
		{
		for (j in (i+1):(length(y))) #Loop on the entries of y	
			{
				mx <- (c(x[i], y[i], x[j], y[j]))	# We have x and y - this is how mx looks in this conf. <---Need to verify !!!!! GW
			####mx <- (c(data[i,m], data[j,m], data[i,n], data[j,n]))	##<-----This is how it was before
				if (length(unique(mx)) >= 2) 							#These are the core calculations 
					{
					if (((mx[1]<mx[3])&(mx[2]<mx[4])) | ((mx[1]>mx[3])&(mx[2]>mx[4]))) 
						{ cosum <- cosum + 1 }
					if (((mx[1]>mx[3])&(mx[3]<mx[4])&(mx[4]>mx[2])&(mx[2]<mx[1])) | ((mx[1]<mx[3])&(mx[3]>mx[4])&(mx[4]<mx[2])&(mx[2]>mx[1]))) 
						{ cesum <- cesum + 1 }  
					}
				#####   Debug===> cat('\n In the comparisons',' i=',i,' j=',j,  'x=',x,' y=',y,' ===>mx=',mx,'cesum and cosum ==',cesum, cosum,'\n') 

			}
		}
	#####   Debug===> cat('\nAfter calculations\n')
	#####   Debug===> cat('\n After comparisons  cesum and cosum ==',cesum, cosum,'\n')
	cesum_adj <- cesum * adj			
	ijsum <- (cosum - cesum_adj)
	return(ijsum)				
}