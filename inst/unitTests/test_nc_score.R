test.nc.score <- function()
  #################################################################################
  #   test.nc.score                                                               #
  #   Function to verify  nc.score                                                #
  #################################################################################
{
    library(RUnit)
    library(infotheo)
	
	data <- read.table("nc_score_input_test.txt",header=TRUE,row.names=1)
	
	
	tol = 0.1
 
	context("NC-score")
	checkEqualsNumeric(nc.score( x=c(3, 1, 1, 3 ,2 ,3, 3 ,1 ,2, 1),  y= c(1 ,2,3 ,3, 1 ,1, 3 ,3, 1, 2)), -0.25,tolerance=tol) 

	testdata<-matrix(c(0.29787234, 0.2978723, 0.2553191, 0.1489362,
		0.17073171, 0.3170732, 0.2682927, 0.2439024,
		0.09302326, 0.3255814, 0.2558140, 0.3255814,
		0.32352941, 0.3235294, 0.1470588, 0.2058824,
		0.17241379, 0.1724138, 0.4137931, 0.2413793,
		0.29729730, 0.2162162, 0.2702703, 0.2162162,
		0.22500000, 0.3250000, 0.2000000, 0.2500000,
		0.12820513, 0.3589744, 0.2307692, 0.2820513,
		0.20000000, 0.2250000, 0.2250000, 0.3500000,
		0.10256410, 0.3076923, 0.1794872, 0.4102564
		),nrow=10,ncol=4,byrow = TRUE)
	dimnames(testdata) = list(
            c("Subject 1", "Subject 2","Subject 3","Subject 4","Subject 5","Subject 6","Subject 7","Subject 8","Subject 9","Subject 10"),
            c("bug 1", "bug 2", "bug 3","bug 4")) # column names 
		
	
	
	
	nc.score.results <-nc.score( x=testdata)
	nc.score.predicted.results <-matrix(c(NA ,-0.25000, -0.21875, -0.65625,
                                       -0.25000 ,      NA, -0.21875 , 0.34375,
                                     -0.21875, -0.21875 ,      NA, -0.21875,
                                     -0.65625,  0.34375, -0.21875,       NA),
                                    nrow=4,ncol=4,byrow = TRUE)
	checkEqualsNumeric(nc.score.predicted.results, nc.score.results, tolerance = tol) 






	
	# Checking two vector version
	mybins <- c(0,0.0001,0.1,0.25)
	x <- data[,1]
	y <- data[,2]
	x.disc <- findInterval(x,mybins)
	y.disc <- findInterval(y,mybins)
	x.disc2 <- discretize(x[which(!is.na(x))])[,1]
	y.disc2 <- discretize(y[which(!is.na(y))])[,1]
	
	
	# Are the values numerically equivalent to Kendall's tau on discretized data?
	checkEqualsNumeric( nc.score(data[,1],data[,2],mybins),cor(x.disc,y.disc,method="kendall",use="complete"), tolerance = tol)
	checkEqualsNumeric( nc.score(data[,1],data[,2]), cor(x.disc2,y.disc2,method="kendall",use="complete"), tolerance = tol)
	
	
	
	# is nc.score.result in the output?
	checkIdentical ( "nc.score.result" %in% names(nc.score(data[,1],data[,2],mybins,verbose=TRUE)),TRUE)
	checkIdentical( "nc.score.result" %in% names(nc.score(data[,1],data[,2],verbose=TRUE)), TRUE)
	
	# Is the number of bins 4?
	checkIdentical( nc.score(data[,1],data[,2],mybins,verbose=TRUE)$n.bins,  4)
	
	# Is the number of bins correctly specified?   
	checkIdentical(nc.score(data[,1],data[,2],verbose=TRUE)$n.bins, floor(sqrt(sum(!is.na(x)))) )
	
	
	# Checking matrix version w/o filtered columns
	####  What is the check here? 
	mybins <- c(0,0.0001,0.1,0.25)
	x <- data
	x.nomiss <- na.omit(data)
	x.disc <- apply(x.nomiss,2,findInterval,vec=mybins)
	x.disc2 <- discretize(x.nomiss)
	nc.score(data[,c(1,2)],bins=mybins,min.samples=0.01,verbose=TRUE)

	# Are the numeric values equal?
	checkEqualsNumeric(0,sum(nc.score(data,bins=mybins,min.samples=0.01)!=cor(x.disc,method="kendall",use="complete"),na.rm=TRUE))
	
	# Is only the diagonal missing?
	checkEqualsNumeric(0,sum(
	  which(
		is.na(nc.score(data,bins=mybins,min.samples=0.01))
	  )
	  !=intersect(
		which(
		upper.tri(
			nc.score(data,bins=mybins,min.samples=0.01),
			diag=TRUE
		)
		),
		which(
		lower.tri(
			nc.score(data,bins=mybins,min.samples=0.01),
			diag=TRUE
		)
		)
	  )
	))  

	# Are numeric values equal?
	checkEqualsNumeric(0,
	sum(
	  nc.score(data,min.samples=0.01) != cor(x.disc2,method="kendall",use="complete"),
	  na.rm=TRUE
	))
	
	# Is only the diagonal missing?
	checkEqualsNumeric(0,
		sum(
		  which(is.na(nc.score(data,min.samples=0.01)))
		  !=intersect(
			which(upper.tri(nc.score(data,min.samples=0.01),diag=TRUE)),
			which(lower.tri(nc.score(data,min.samples=0.01),diag=TRUE))
		  )
		))
	
	
	## Checking matrix version w/ filtered columns
	filter_cols <- which(
	  apply(
		data,
		2,
		function(col){
		length(which(col > 0.0001))/sum(!is.na(col))
		}
	  ) < 0.2
	)

	
	
			
			
		# Are numeric values equal?
		z <- nc.score(data,verbose=TRUE,min.samples=0.2)
		checkEqualsNumeric(0,
			sum(
			  z$nc.score.matrix!=cor(x.disc2,method="kendall",use="complete"),
			  na.rm=TRUE
			)
		)
		
		# Are the missing values correct?
		checkEqualsNumeric(0,
			sum(
			  which(is.na(z$nc.score.matrix))
			  !=sort(
				Reduce(
				union,
				list(
					intersect(
					which(upper.tri(z$nc.score.matrix,diag=TRUE)),
					which(lower.tri(z$nc.score.matrix,diag=TRUE))
					),
					sapply(nrow(z$nc.score.matrix)*(z$columns.not.passing.qc-1),'+',(1:nrow(z$nc.score.matrix))),
					sapply(z$columns.not.passing.qc,seq,to=ncol(z$nc.score.matrix)*nrow(z$nc.score.matrix),by=nrow(z$nc.score.matrix))
				)
				)
			  )
			)
			)
		# Should give an error
		checkException(nc.score(data,y))

		
		

		# Should give an error
		checkException(nc.score(data[,1],data[1:20,2]))

		# Should give an error
		checkException(nc.score(bins=mybins))

		# Should give an error
		checkException(nc.score(y))
		
		
		

		# Should give a warning about replacing value - This is just a warning 
		options(warn=2)               #Set warning level to 2 so warnings become errors
		checkException(nc.score(data[,1],data[,2],bins=-1))

		# Should give a warning about replacing value - This is just a warning 
		checkException(nc.score(data[,1],data[,2],bins="4"))

		# Should give a warning about replacing value - This is just a warning 
		checkException(nc.score(data[,1],data[,2],bins=pi))

		# Should give a warning about replacing value
		checkException(nc.score(data,min.samples=1.5))

		# Should give a warning about replacing value
		checkException(nc.score(data,min.samples=-2))

		# Should give a warning about replacing value
		checkException(nc.score(data,min.samples="0.1"))

		# Should give a warning about replacing value
		checkException(nc.score(data,min.samples=c(0.1,0.2)))

		
		# Should give a warning about replacing value
		checkException(nc.score(data,min.abundance="0.1"))

		# Should give a warning about replacing value
		checkException(nc.score(data,min.abundance=c(0.0001,0.0002)))
	
	
		# Should give a warning about replacing value
		#** No warning, value returned
		checkException(nc.score(data[,1],data[,2],min.samples=1.5))

		# Should give a warning about replacing value
		#** No warning, value returned
		checkException(nc.score(data[,1],data[,2],min.samples=-2))

		# Should give a warning about replacing value
		#** No warning, value returned
		checkException(nc.score(data[,1],data[,2],min.samples="0.1"))

		# Should give a warning about replacing value
		#** No warning, value returned
		checkException(nc.score(data[,1],data[,2],min.samples=c(0.1,0.2)))

		# Should give a warning about replacing value
		#** No warning, value returned
		checkException(nc.score(data[,1],data[,2],min.abundance=c(0.0001,0.0002)))

	
		options(warn=0)   #Reset  warning level to 0
  }
 