\name{nc.score}
\alias{nc.score}
\title{
  nc.score
}
\description{
  Calculates species-level co-variation and co-exclusion patterns based on an extension of the checkerboard score to ordinal data.
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
   First entity to be processed.
}
  \item{y}{
   Second entity to be processed.
}
  \item{bins}{
 Number of bins to use
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
\examples{
	data1<-matrix(rlnorm(21*5,1,1),nrow=21)
	nc.score.results <- nc.score( x=data1 )
}


