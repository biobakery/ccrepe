\name{nc.score}
\alias{nc.score}
\title{
    Calculates species-level co-variation and co-exclusion patterns based on an extension of the checkerboard score to ordinal data.
}
\description{
  NC score ordinalizes relative abundance data by "binning."  The resulting ordinal data is evaluated for the presence of cooccurrence and coexclusion patterns.
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
  Minimum abundance threshold: Bugs which do not have at least min.samples samples with abundance min.abundance are excluded to reduce computational burden.
}
  \item{min.samples}{
  Minimum samples threshold: Bugs which do not have at least min.samples samples with abundance min.abundance are excluded to reduce computational burden.
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


