test.nc.score <- function()
  #################################################################################
  #   test.nc.score                                                               #
  #   Function to verify  nc.score                                                #
  #################################################################################
{
	tol = 0.1
	
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

	z <- nc.score(data,bins=mybins,min.samples=0.2,verbose=TRUE)

	# Correct columns filtered?
	checkEqualsNumeric(z$columns.not.passing.qc,  filter_cols)

	
	
	
	
  }
 