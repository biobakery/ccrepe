nc.score.helper <- function(data,CA) {
	#****************************************************************************************
	#*   Core calculations                                                                  *
	#****************************************************************************************
	mode(data) <- "numeric"
	CA$nc.score.matrix <- matrix(nrow=ncol(data),ncol=ncol(data))		#Build results area
	n <- length(unique(c(data)))
	adj <- ((1.5)*n*(n-1)/(n^2-n+1))
	for(i in 1:(ncol(data))) {    
		for(j in (i):ncol(data)) {   		
			ijsum <- 0
			cosum <- 0
			cesum <- 0
			for (m in 1:(nrow(data))) {   
				for (n in (m):nrow(data)) {  	  
					mx <- (c(data[m,i], data[m,j], data[n,i], data[n,i]))
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
		CA$nc.score.matrix[i,j] <- ijsum				#Post it to the matrix
		CA$nc.score.matrix[j,i] <- ijsum				#Post it to the matrix
		}
	}
  return(CA$nc.score.matrix)							#Return the calculated results matrix
}
 