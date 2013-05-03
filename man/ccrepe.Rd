\name{ccrepe}
\alias{ccrepe}
\title{
  ccrepe 
   A function to calculate correlation, p-values and q-values of the relative abundance of bugs in one or two cavities
}
\description{
 
   The function uploads the relative abundances of bugs in one or two cavities and calculates their correlation, p-value and q-value ussing permutation and 
   rebooting of the data
}
\usage{
 ccrepe(x = NA, 
 y = NA, 
 method = 'pearson',
 min.subj = 20,  
 iterations = 1000,
 subset.cols.1 = c(0),
 subset.cols.2 = c(0),
 errthresh = 0.0001,
 distributions = NA)
}
\arguments{
  \item{x}{
   First dataframe containing the relative abundances in cavity1 : Columns are bugs, rows are samples
}
  \item{y}{
   Second dataframe containing the relative abundances in cavity2: Columns are bugs, rows are samples
}
 
  \item{method}{
	Type of correlation: Can be one of three: pearson, spearman or kendall (default: pearson)
}
  \item{min.subj}{
  Minimum rows (Subjects) in "data" frame to proceed (default: 20)
}
  \item{iterations}{
  Reboot iterations (default: 1000)
}
  \item{subset.cols.1}{
  Subset of cols(bugs) from cavity1 to iterate on (default: c(0)== ALL)
}
  \item{subset.cols.2}{
  Subset of cols(bugs) from cavity2 to iterate on (default: c(0)== ALL)
}
  \item{errthresh}{
	Threshold error if there is enough data to calculate cor an pval for certain i and k (default: 0.0001)
}
  \item{distributions}{
	Output Distribution file  (default:NA)
}
}
 
\value{

	Returns a dataframe containing the calculations results and the parameters used.\cr
	Parameters shown:\cr
		min.subj, iterations, errthresh, method, subset.cols.1  (and subset.cols.2 -if applicable)\cr
		Their descriptions are above\cr
	The calculation results consist of:\cr
		PValues.matrix,  Q.matrix and  data.cor that contains:  bug1 bug2 cor   p.value  q.value
}
\author{
  Emma Schwager <emma.schwager@gmail.com>
}
\examples{
	mydata<-matrix(rlnorm(21*15,1,1),nrow=21)
	ccrepe.results <- ccrepe (x=mydata,iterations=100)
}
