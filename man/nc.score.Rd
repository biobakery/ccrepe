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
 Number of bins to use
 The number of bins are used to discretize the input matrix.
 The package uses the discretize function of the package infotheo where the number of bins is an integer specifying the number of bins to be used for the discretization 
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


  Diamond JM (1975) Assembly of species communities. In: Cody ML, Diamond JM,
editors. Ecology and Evolution of Communities. Cambridge, MA: Harvard
University Press. pp. 342-444.

}
\examples{
	data1<-matrix(rlnorm(21*5,1,1),nrow=21)
	nc.score.results <- nc.score( x=data1 )
}


