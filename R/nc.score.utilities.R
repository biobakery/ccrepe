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

}

process_wrapper_input <- function(CA,
                                     input.min.abundance,
                                     input.min.samples){
        CA <-list()		
	CA <- process_min_abundance_min_samples( CA,                  #Process minimum abundance and min samples
			input.min.abundance,
			input.min.samples)								 

 
	#*********************************************************
	#* Filter                                                *
	#*********************************************************
        CA <- qc_filter('x',CA)
        CA <- qc_filter('y',CA)

	return(CA)

}

qc_filter <- function(var_name,CA){
    var_values <- CA[[var_name]]
    if(is.vector(var_values)){
        non.miss <- var_values[-which(is.na(var_values))]
        if( sum( non.miss >= CA$min.abundance ) <= CA$min.samples*length(non.miss) ){
            stop(paste0(
                "'",var_name,"' has ",
                100*(sum(non.miss >= CA$min.abundance)/length(non.miss)),
                "% of non-missing samples with abundance less than min.abundance=",
                CA$min.abundance
                ))
        } else {
            CA[[paste0(var_name,".filtered")]] <- var_values
            CA[[paste0(var_name,".features.filtered")]] <- c()
        }
    } else if (is.matrix(var_values) || is.data.frame(var_values)){
        CA[[paste0(var_name,".features.filtered")]] <- which(
            apply(x,2,
                  function(col){
                      sum(col[which(!is.na(col))] >= CA$min.abundance)/length(col[which(!is.na(col))])
                  }
                  ) <= CA$min.samples
            )
        if(length(CA[[paste0(var_name,".features.filtered")]])==ncol(var_values)){
            stop(paste0(
                "'",var_name,"' has 0 features which have 100*min.samples=",
                100*CA$min.samples,"% of samples with abundances >=",
                CA$min.abundance, "(min.abundance)."
                ))
        }
        CA[[paste0(var_name,".filtered")]] <- var_values[,-CA[[paste0(var_name,".features.filtered")]]]
    } else if (is.null(var_values)){ 
        CA[[paste0(var_name,".filtered")]] <- var_values
        CA[[paste0(var_name,".features.filtered")]] <- c()
    } else {
        stop(paste0("Unrecognized input type for variable '",var_name,"'"))
    }
    return(CA)
    
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


