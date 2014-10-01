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
                .Call(C_cor,
                      Bin(dropNA(x, nas),nbins,bin.cutoffs),
                      Bin(dropNA(y, nas),nbins,bin.cutoffs),
                      na.method,
                      method == "kendall")
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


## nc.score.wrapper <-
## function(
## #*************************************************************************************
## #*  	nc.score                                                                     *
## #*************************************************************************************
##    	x,						#First input 
## 	y=NULL,						#Second input
##         use = "everything",                             # What to do with missing values
## 	nbins=NULL,					#Number of Input Bins
##         bin.cutoffs=NULL,                               #Cutoff values for the bins
## 	verbose = FALSE,			#Request for verbose output?
## 	min.abundance=0.0001,		#Minimum Abundance
## 	min.samples=0.1)			#Minimum Samples
	
## {
## #*************************************************************************************
## #*  	NC_Score                                                                     *
## #*************************************************************************************

##     na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs",
##                                "everything", "na.or.complete"))

##     processed_input <- preprocess_nc_score_input(
##         x,
##         y,
##         na.method,
##         nbins,
##         bin.cutoffs
##         )

##     x <- processed_input$x
##     y <- processed_input$y

##     nbins       <- processed_input$nbins
##     bin.cutoffs <- processed_input$bin.cutoffs

##     CA <- list()									#Set the common area
##     CA$x <- x
##     CA$y <- y
##     CA$use <- use
##     CA$nbins <- nbins
##     CA$bin.cutoffs <- bin.cutoffs
##     CA$verbose     <- verbose
    
##     CA <- process_wraper_input(CA,min.abundance,min.samples)

##     CA$nc.score.result <- nc.score(
##         CA$x.filtered,
##         CA$y.filtered,
##         CA$use,
##         CA$nbins,
##         CA$bin.cutoffs
##         )

##     if(CA$verbose){
##          return(CA)
##     } else {
##         return(CA$nc.score.result)
##     }
## }

