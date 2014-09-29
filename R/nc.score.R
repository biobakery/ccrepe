nc.score <- function (x, y = NULL, use = "everything",nbins=NULL,bin.cutoffs=NULL)
    {
        method = "kendall"
        C_cor = get("C_cor",asNamespace("stats"))
        na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs",
                                   "everything", "na.or.complete"))

        processed_input <- preprocess_nc_score_input(
            x,
            y,
            na.method,
            nbins,
            bin.cutoffs
            )

        x <- processed_input$x
        y <- processed_input$y
        
        nbins       <- processed_input$nbins
        bin.cutoffs <- processed_input$bin.cutoffs

        if (na.method %in% c(2L, 5L)) {
            if (is.null(y)) {
                .Call(C_cor, Bin(na.omit(x),nbins,bin.cutoffs), NULL, na.method, method ==
                      "kendall")
            }
            else {
                nas <- attr(na.omit(cbind(x, y)), "na.action")
                dropNA <- function(x, nas) {
                    if (length(nas)) {
                        if (is.matrix(x))
                            x[-nas, , drop = FALSE]
                        else x[-nas]
                    }
                    else x
                }
                .Call(C_cor, Bin(dropNA(x, nas),nbins,bin.cutoffs), Bin(dropNA(y,
                                                                               nas),nbins,bin.cutoffs), na.method, method == "kendall")
            }
        }
        else if (na.method != 3L) {
            x <- Bin(x,nbins,bin.cutoffs)
            if (!is.null(y))
                y <- Bin(y,nbins,bin.cutoffs)
            .Call(C_cor, x, y, na.method, method == "kendall")
        }
        else {
            if (is.null(y)) {
                ncy <- ncx <- ncol(x)
                if (ncx == 0)
                    stop("'x' is empty")
                r <- matrix(0, nrow = ncx, ncol = ncy)
                for (i in seq_len(ncx)) {
                    for (j in seq_len(i)) {
                        x2 <- x[, i]
                        y2 <- x[, j]
                        ok <- complete.cases(x2, y2)
                        x2 <- Bin(x2[ok],nbins,bin.cutoffs)
                        y2 <- Bin(y2[ok],nbins,bin.cutoffs)
                        r[i, j] <- if (any(ok))
                            .Call(C_cor, x2, y2, 1L, method == "kendall")
                        else NA
                    }
                }
                r <- r + t(r) - diag(diag(r))
                rownames(r) <- colnames(x)
                colnames(r) <- colnames(x)
                r
            }
            else {
                if (length(x) == 0L || length(y) == 0L)
                    stop("both 'x' and 'y' must be non-empty")
                matrix_result <- is.matrix(x) || is.matrix(y)
                if (!is.matrix(x))
                    x <- matrix(x, ncol = 1L)
                if (!is.matrix(y))
                    y <- matrix(y, ncol = 1L)
                ncx <- ncol(x)
                ncy <- ncol(y)
                r <- matrix(0, nrow = ncx, ncol = ncy)
                for (i in seq_len(ncx)) {
                    for (j in seq_len(ncy)) {
                        x2 <- x[, i]
                        y2 <- y[, j]
                        ok <- complete.cases(x2, y2)
                        x2 <- Bin(x2[ok],nbins,bin.cutoffs)
                        y2 <- Bin(y2[ok],nbins,bin.cutoffs)
                        r[i, j] <- ifelse(
                            (any(ok)),
                            .Call(C_cor, x2, y2, 1L, method == "kendall"),
                            NA )
                    }
                }
                rownames(r) <- colnames(x)
                colnames(r) <- colnames(y)
                if (matrix_result){
                    r
                } else drop(r)
            }
        }
    }


## nc.score <-
## function(
## #*************************************************************************************
## #*  	nc.score                                                                     *
## #*************************************************************************************
##    	x=NULL,						#First input 
## 	y=NULL,						#Second input
## 	bins=NA,					#Number of Input Bins
## 	verbose = FALSE,			#Request for verbose output?
## 	min.abundance=0.0001,		#Minimum Abundance
## 	min.samples=0.1)			#Minimum Samples
	
## {
## #*************************************************************************************
## #*  	NC_Score                                                                     *
## #*************************************************************************************
 
## CA <- list()									#Set the common area
## if(is.data.frame(x)) x <- as.matrix(x)
## if (!is.matrix(x) && is.null(y)){
##     stop("supply both 'x' and 'y' or a matrix-like 'x'")
## }
## if (!(is.numeric(x) || is.logical(x))){
##     stop("'x' must be numeric")
## }
## if(!is.null(y)){
##     if(is.matrix(x)) stop("supply matrix-like 'x' only if 'y' not supplied")
##     if(length(x) != length(y)) stop('incompatible dimensions')
##     nonmissing.subjects <- intersect(which(!is.na(x)),which(!is.na(y)))      #Find subjects which are present in both x and y
##     CA$x <- x[nonmissing.subjects]					#Post x without missing subjects to common area
##     CA$y <- y[nonmissing.subjects]					#Post y without missing subjects to common area
##     CA$verbose <- verbose						#Post the verbose flag


##     CA <- process.input.bins(bins, CA)			#Process bins input request
 
## 	CA <- process_min_abundance_min_samples( CA,                  #Process minimum abundance and min samples
## 			min.abundance,
## 			min.samples)
	
	
	
## 	CA$xy <- data.frame(CA$x, CA$y)             #Merge the two columns into a data frame for filtering
	
	
## 	CA <- qc_filter(CA$xy ,CA)
  
	
##     if (is.null(CA$bins))					#Check if bins is a number or a vector with entries
##         {
##             CA$x.discretized = discretize(CA$x,nbins = CA$n.bins)[,1]	#Discretize x
##             CA$y.discretized = discretize(CA$y,nbins = CA$n.bins)[,1]	#Discretize y
##         } else {
##             CA$x.discretized  = findInterval(CA$x,CA$bins)	 #Use the bins provided by the User
##             CA$y.discretized  = findInterval(CA$y,CA$bins)	 #Use the bins provided by the User
##         }
 
	
	
##     CA$nc.score.result <- cor(CA$x.discretized,CA$y.discretized,method="kendall")
##     if (CA$verbose == TRUE)
##         {return(CA)}
##     else
##         {return(CA$nc.score.result)}
## } else { 
		
## 	#************************************************************************
## 	#*        It is a dataframe or a matrix                                 *
## 	#************************************************************************

##     CA = preprocess_nc_score_input (
##         x, 										#First Input
##         bins,									#Bins
##         min.abundance,							#Minimum Abundance
##         min.samples)							#Minimum Samples

##     CA$verbose <- verbose							#Post the verbose flag


##     if (is.null(CA$bins))					#Check if bins is a number or a vector with entries
##         {
##             CA$x.discretized = discretize(CA$x.filtered,nbins=CA$n.bins)	#Bins is a number; Discretize it and post the value in the discretized matrix
##         }
##     else
##         {
##             CA$x.discretized = apply(CA$x.filtered,2,findInterval,vec=CA$bins)	 #Use the bins provided by the User
##         }


##     CA$nc.score.matrix <- cor(CA$x.discretized,method="kendall")

##     diag(CA$nc.score.matrix)<-NA	#We are setting the diagonal entries in the matrix to NA
##     if (length(CA$columns.not.passing.qc) > 0)  #If there were columns that did not pass QA, we need to add corr cols with NA
##         {
##             original.nc.score.dim <- ncol(CA$nc.score.matrix)  #Columns in the original matrix
##             rebuilt.matrix <- CA$nc.score.matrix					#Allocate the rebuilt matrix


##             for (indx in 1:length(CA$columns.not.passing.qc))
##                 {

##                     if (CA$columns.not.passing.qc [indx] -1 <  ncol(rebuilt.matrix) )
##                         {
##                             rebuilt.matrix <- cbind(rebuilt.matrix[,1:CA$columns.not.passing.qc [indx] -1],		#Left part of the rebuilt matrix
##                                                     rep(NA, original.nc.score.dim),   #Insert NA's
##                                                     rebuilt.matrix[,CA$columns.not.passing.qc [indx]:ncol(rebuilt.matrix)])		#Right part of the rebuilt matrix
##                         }
##                     else
##                         {
##                             rebuilt.matrix <- cbind(rebuilt.matrix[,1:CA$columns.not.passing.qc [indx] -1],		#Left part of the rebuilt matrix
##                                                     rep(NA, original.nc.score.dim))		#Insert NA's
##                         }
##                 }

##             for (indx in 1:length(CA$columns.not.passing.qc))
##                 {
##                     if (CA$columns.not.passing.qc [indx] -1 <  nrow(rebuilt.matrix) )
##                         {
##                             rebuilt.matrix <- rbind(rebuilt.matrix[1:CA$columns.not.passing.qc [indx] -1,],				#Upper part of the rebuilt matrix
##                                                     rep(NA, ncol(rebuilt.matrix)),   #Insert NA's
##                                                     rebuilt.matrix[CA$columns.not.passing.qc [indx]:nrow(rebuilt.matrix),])		 	##Lower Part of the rebuilt matrix
##                         }
##                     else
##                         {
##                             rebuilt.matrix <- rbind(rebuilt.matrix[1:CA$columns.not.passing.qc [indx] -1,],				#Upper part of the rebuilt matrix
##                                                     rep(NA, ncol(rebuilt.matrix)))   #Insert NA's
##                         }
##                 }

##             colnames(rebuilt.matrix)<- CA$original.column.names   #Post the original column names into the column names of rebuilt matrix
##             rownames(rebuilt.matrix)<- CA$original.column.names   #Post the original column names into the row names
##             CA$nc.score.matrix <- rebuilt.matrix				#Post the matrix
##         }

## #    CA$x.discretized <- NULL		#Not needed anymore
## #    CA$x <- NULL					#Not Needed anymore
##     CA$x.filtered <- NULL
##     CA$input.total.cols <- NULL		#Not needed anymore
##     if(length(CA$columns.not.passing.qc) == 0){
##         CA$columns.not.passing.qc <- NULL  #Not needed anymore
##     }
##     CA$original.column.names <- NULL  #Not needed anymore
##     CA$names.of.cols.failing.qc <- NULL  # Not needed anymore

##     if (!CA$verbose == TRUE)		#If abbreviated output
##         {
##             CA <- CA$nc.score.matrix	#just post the resulting matrix
##         }
## }
    
## return(CA)
## }
