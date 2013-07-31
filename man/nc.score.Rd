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
   First entity to be processed.   Columns are bugs, rows are samples.
 }
  \item{y}{
   Vector.  If selected, x must be a vector as well.   
}
  \item{bins}{
 Either the number of bins to use or a vector specifying the cutoffs.
 If a single number is give, this is used as the number of bins in the discretize function of the package infotheo.
 If a vector is specified, the function findInterval is used to discretize the data according to the cutoffs given.
 The default behavior is to use the defaults for the discretize function.
}
\item{verbose}{
 Request verbose output
}
  \item{min.abundance}{
  Minimum abundance threshold: Assures selection of species with abundance >= min.abundance in more than min.samples percent of samples.
}
  \item{min.samples}{
  Minimum samples threshold: Assures selection of species with abundance >= min.abundance in more than min.samples percent of samples.
}
}
\value{ 
	Matrix or vectors of normalized scores.
}

\author{
 Craig Bielski<craig.bielski@gmail.com>
}
\references{


  First Author et al (2014) Detecting statistically significant associtations between sparse and high dimensional compositioanl data.

}
\examples{	
        data1<-matrix(c(
        0.29787234, 0.2978723, 0.2553191, 0.1489362,
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

	nc.score.results <- nc.score( x=data1 )
        nc.sccore.results.bins <- nc.score( x=data1, bins=5 )
        nc.score.results.bin.cutoffs <- nc.score( x=data1, bins = c(.001,.25,.5) )
}


