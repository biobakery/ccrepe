\name{nc.score}
\alias{nc.score}
\title{
  Calculate  species-level co-variation and co-exclusion patterns based on an extension of the checkerboard score to ordinal data

}
\description{
  Calculates species-level co-variation and co-exclusion patterns based on an extension of the checkerboard score to ordinal data
  }
\usage{
nc.score(x = NA,
 y = NA, 
 adBins = NA, 
 min.abundance = 1e-04,
 min.samples = 0.1)
}
\arguments{
  \item{x}{
   First data frame to be processed
}
  \item{y}{
   Vector.  
   If x is entered then both x and t must be only vectors and the program calculates the score as a number
}
  \item{adBins}{
 Number of bins to use
}
  \item{min.abundance}{
  Minimum abundance threshold: Assures selection of species with abundance >= min.abundance in more than min.samples percent of samples
}
  \item{min.samples}{
  Minimum samples threshold: Assures selection of species with abundance >= min.abundance in more than min.samples percent of samples
}
}
\value{ 
	matrix or vectors of normalized scores
}

\author{
 Craig Bielski<craig.bielski@gmail.com>
}
\examples{

	x1<-matrix(rlnorm(10*9,1,1),nrow=10)
	y1<-matrix(rlnorm(10*9,1,1),nrow=10)
	nc.score.result <-nc.score (x=x1,y=y1)
}


