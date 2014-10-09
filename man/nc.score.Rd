\name{nc.score}
\alias{nc.score}
\title{


  nc.score
  
  
}
\description{


  nc.score calculates species-level co-variation and co-exclusion patterns based on an extension of the checkerboard score to ordinal data.
  
  It is an extension to Diamond's checkerboard score (See references below) to ordinal data and implements a
  framework for robust detection of species-level association patterns in metagenomic data.
  It performs quality control and transforms relative abundances to ordinal values and computes an association measure,
  the NC-score (n-dimensional checkerboard score), for each species pair.
  
  }
\usage{
nc.score(x,
 y = NULL, 
 use = "everything",
 nbins = NULL, 
 bin.cutoffs = NULL
)
}
\arguments{
  \item{x}{
   A numeric vector, data frame, or matrix.  The first entity to be processed.   Columns are bugs, rows are samples.
 }
  \item{y}{
   Numeric vector.  If selected, x must be a numeric vector as well.   
}
	\item{use}{
   use = "everything"
	}
	\item{nbins}{
   The number of bins
	}
	\item{bin.cutoffs}{
   bin.cuttofs 
	}
}
\value{ 
	Matrix or vector of normalized scores.
}

\author{
 Craig Bielski<craig.bielski@gmail.com>
}
\references{


  Emma Schwager and Colleagues.  Detecting statistically significant
  associtations between sparse and high dimensional compositioanl
  data. In Progress.

}
\examples{	
data <- matrix(rlnorm(40,meanlog=0,sdlog=1),nrow=10)
data.rowsum <- apply(data,1,sum)
data.norm <- data/data.rowsum
testdata <- data.norm
dimnames(testdata) <- list(paste("Sample",seq(1,10)),paste("Feature",seq(1,4)))

nc.score.results <- nc.score( x=testdata )
nc.score.results.bins <- nc.score( x=testdata  )
nc.score.results.bin.cutoffs <- nc.score( x=testdata  )
nc.score.results
nc.score.results.bins
nc.score.results.bin.cutoffs
}


