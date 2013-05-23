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
 bins = NA, 
 min.abundance = 1e-04,
 min.samples = 0.1)
}
\arguments{
  \item{x}{
   First entity to be processed
}
  \item{y}{
   Second entity to be processed
}
  \item{bins}{
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
	data1<-matrix(rlnorm(21*5,1,1),nrow=21)
	data2<-matrix(rlnorm(5*17,1,1),nrow=17)
	nc.score.results <- ccrepe( x=data1, y=data2, method = nc.score )
}


