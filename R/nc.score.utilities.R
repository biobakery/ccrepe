Bin <- function(u,nbins,bin.cutoffs) {
    if(!is.null(nbins)){
        if (length(u) == 0L)
            u
        else if (is.matrix(u)) {
            if (nrow(u) > 1L)
                as.matrix(discretize(u,nbins=nbins))
            else as.matrix(discretize(u,nbins=nbins))
        }
        else discretize(u,nbins=nbins)[,1]
    } else if (!is.null(bin.cutoffs)){
        if (length(u) == 0L)
            u
        else if (is.matrix(u)){
            if(nrow(u) > 1L)
                apply(u,2,findInterval,vec=bin.cutoffs)
            else matrix( apply(u,2,findInterval,vec=bin.cutoffs),nrow=1 )
        }
        else findInterval(u,vec=bin.cutoffs)
    }
}

preprocess_nc_score_input <-
function(
#********************************************************************************************
#*  	Pre process input and build the common area                                         *
#********************************************************************************************
    x,
    y,
    na.method,
    nbins,
    bin.cutoffs
    )
{
      if (is.na(na.method))
          stop("invalid 'use' argument")
      if (is.data.frame(y))
          y <- as.matrix(y)
      if (is.data.frame(x))
          x <- as.matrix(x)
      if (!is.matrix(x) && is.null(y))
          stop("supply both 'x' and 'y' or a matrix-like 'x'")
      if (!(is.numeric(x) || is.logical(x)))
          stop("'x' must be numeric")
      stopifnot(is.atomic(x))
      if (!is.null(y)) {
          if (!(is.numeric(y) || is.logical(y)))
              stop("'y' must be numeric")
          stopifnot(is.atomic(y))
      }
      if(!is.null(nbins) && !is.null(bin.cutoffs)){
          stop("supply only one of 'nbins' or 'bin.cutoffs'")
      }
      if(!is.null(nbins) && !is.numeric(nbins)){
          stop("supply a numeric bin number")
      }
      if(!is.null(nbins) && length(nbins) > 1L){
          warning("nbins has length > 1 and only the first element will be used")
          nbins <- nbins[1]
      }
      if(!is.null(nbins) && nbins <= 0){
          stop("supply a positive bin number")
      }
      if(!is.null(nbins) && (nbins %% 1 != 0)){
          stop("supply an integer bin number")
      }
      if(!is.null(bin.cutoffs) && !is.numeric(bin.cutoffs)){
          stop("invalid bin cutoff values")
      }
      if(!is.null(bin.cutoffs) && sum(sort(bin.cutoffs)!=bin.cutoffs)>0){
          warning("bin.cutoffs not monotonically ordered - using sorted values")
          bin.cutoffs <- sort(bin.cutoffs)
      }
      if(is.null(nbins) && is.null(bin.cutoffs)){
          if(is.matrix(x)){
              nbins <- floor(sqrt(nrow(x)))
          } else {
              nbins <- floor(sqrt(length(x)))
          }
      }
      return(list(
          x=x,
          y=y,
          nbins=nbins,
          bin.cutoffs=bin.cutoffs
          ))

        ## CA <-list()		
	## CA <- process_min_abundance_min_samples( CA,                  #Process minimum abundance and min samples
	## 		input.min.abundance,
	## 		input.min.samples)								 

 
	## #*********************************************************
	## #* Filter                                                *
	## #*********************************************************
	## CA$x <- x												#Post x in the common area
        ## CA$x <-na.omit(CA$x)											#remove NAs
	## CA <-qc_filter(CA$x,CA)										#filter by abundance thresholds
	## if (is.null(rownames(CA$x))) {rownames(CA$x)<-seq(1:nrow(CA$x))} #If there are no row names - plug them in
	## #****************************************************************************************
	## #*   Process Bins  passed by the user  or set default                                   *
	## #****************************************************************************************
	## CA <- process.input.bins(input.bins, CA)			#Process Input number of bins or set default
	
	## return(CA)
}



process.input.bins <-
function(input.bins, CA)
#********************************************************************************************
#*  	Preprocess input bins as passed by the user or use default                          *
#********************************************************************************************
{
	if (length(input.bins) > 1)					#If the user passed a set of bins,  sort them and use them
		{
			CA$bins <- sort(input.bins)
			CA$n.bins <- length(intersect(which(CA$bins<1),which(CA$bins>0))) + 1
			if (! 0 %in% CA$bins )				#Check if 0 is in the list of the bins the User passed - not there - add it and sort 
				{ 
				 CA$bins<-c(CA$bins,0)			#Add 0
				 CA$bins <- sort(CA$bins)		#And sort it
				}
			return(CA)
		}
		
	if (class(CA$x) == "numeric") 
		{
		n.bins <- floor(sqrt(length(CA$x)))					#If a vector - take sqrt of the length	
		} else { 
		n.bins <- floor(sqrt(nrow(CA$x))) 					#Set the default number of bins 
		}
		
	if (!is.na(input.bins))								#If User entered number of Bins
		if (!is.numeric(input.bins)){
			warning('Number of bins must be numeric - Using default bins = sqrt(Number of rows) ')
                    } else if (input.bins%%1!=0) {                                       # Check that the number of bins is an integer
			warning('Number of bins must be an integer - Using default bins = sqrt(Number of rows) ')
                    } else if (input.bins <=0){
			warning('Number of bins must be positve - Using default bins = sqrt(Number of rows) ')
                    } else {
			bins = input.bins					#Valid Input from the User - Use it
                    }
	
	CA$bins = NULL										#Post it to Common Area
	CA$n.bins = n.bins									#Claculate the number of bins
	return(CA)
}
qc_filter <-
function(data,CA) {
#********************************************************************************************
#*  	quality control: select species with >= 1E-4 relative abundance in >= 10% of samples*
#********************************************************************************************

	tmp <- {}
	names <- {}
	CA$names.of.cols.failing.qc <- {} 			#Names of the cols failing QC
	CA$original.column.names <- colnames(data)	#These are the original column names
	
	CA$input.total.cols <- ncol(data)			#number of cols in the input
	CA$columns.not.passing.qc = vector()		#define a vector to contain the seq number of cols that did not pass qc

	

	for (i in seq_len(ncol(data))) 
	{
		if (length(which(data[,i] >= CA$min.abundance)) >= CA$min.samples*nrow(data))
		{
			tmp <- cbind(tmp, data[,i])
			names <- c(names, colnames(data)[i])
		} else 
		{
			CA$names.of.cols.failing.qc<-c(CA$names.of.cols.failing.qc, colnames(data)[i])
			CA$columns.not.passing.qc <- c(CA$columns.not.passing.qc,i)		#Add the number of the col that did not pass
			warning_msg <- paste("Column ",i,":  ",colnames(data)[i]," - doesn't pass quality control: try adjusting min.samples or min.abundance")
			warning(warning_msg)
		}
	}
  colnames(tmp) <- names					#Post the column names
  rownames(tmp) <- rownames(data)			#Post the data row names
  CA$x.filtered <- tmp								#Post in the common area

  return(CA)								#Return the common Area
}





process_min_abundance_min_samples <-function(   CA,
			input.min.abundance,
			input.min.samples)
#********************************************************************************************
#*  	Process minimum abundance and input samples                                         *
#********************************************************************************************
{
 	if (!is.numeric(input.min.abundance))	#check if the User entered a valid threshold1
		{
		warning('\nMinimum abundance must be numeric - using default=0.0001\n')
		CA$min.abundance = 0.0001
		 }				#If it is not valid - force default
        else if(length(input.min.abundance) > 1)
            {
                warning("More than one minimum abundance value given - using first one")
                CA$min.abundance = input.min.abundance[1]
            }
	else
		{CA$min.abundance = input.min.abundance}	#Else - use it

	if (!is.numeric(input.min.samples))	#check if the User entered a valid threshold1
		{
		warning('\nMinimum samples must be numeric - using default=0.1\n')
		CA$min.samples = 0.1 				#If it is not valid - force default
		}
	else if (input.min.samples <= 0 || input.min.samples >= 1)	#check if the User entered a valid threshold1
		{
		warning('\nMinimum samples must be between 0 and 1 inclusive - using default=0.1\n')
		CA$min.samples = 0.1 				#If it is not valid - force default
		}
        else if (length(input.min.samples) > 1)
            {
                warning("More than one min.samples value given - using first value")
                CA$min.samples <- input.min.samples[1]
            }
	else
		{CA$min.samples = input.min.samples}	#Else - use it

	return(CA)								#Return the common Area
}


