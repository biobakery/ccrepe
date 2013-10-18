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
nc.score(x = NA,
 y = NA, 
 bins = NA, 
 verbose = FALSE,			 
 min.abundance = 1e-04,
 min.samples = 0.1)
}
\arguments{
  \item{x}{
   A numeric vector, data frame, or matrix.  The first entity to be processed.   Columns are bugs, rows are samples.
 }
  \item{y}{
   Numeric vector.  If selected, x must be a numeric vector as well.   
}
  \item{bins}{
 Either a single integer specifying the number of bins to use or a numeric vector specifying the cutoffs.
 If a single number is given, this is used as the number of bins in the discretize function of the package infotheo.
 If a vector is specified, the function findInterval is used to discretize the data according to the cutoffs given.
 The default behavior is to use the defaults for the discretize function.
}
\item{verbose}{
 Logical flag to request verbose output
}
  \item{min.abundance}{
  A numeric value specifying the minimum abundance threshold: Assures selection of species with abundance >= min.abundance in more than min.samples percent of samples.
}
  \item{min.samples}{
  An fraction specifying the minimum samples threshold: Assures selection of species with abundance >= min.abundance in more than min.samples percent of samples. Default=0.1
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
nc.score.results.bins <- nc.score( x=testdata, bins=5 )
nc.score.results.bin.cutoffs <- nc.score( x=testdata, bins = c(.001,.25,.5) )
nc.score.results
nc.score.results.bins
nc.score.results.bin.cutoffs
}


