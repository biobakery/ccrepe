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
ccrepe(frmOne = NA, frmTwo = NA,  method = 'pearson',  min.rows = 20, Iterations = 1000, subset1 = c(0), subset2 = c(0), errthresh = 0.0001, distributions = NA)
}
\arguments{
  \item{frmOne}{
   First dataframe containing the relative abundances in cavity1 : Columns are bugs, rows are samples
}
  \item{frmTwo}{
   Second dataframe containing the relative abundances in cavity2: Columns are bugs, rows are samples
}
 
  \item{method}{
	Type of correlation: Can be one of three: pearson, spearman or kendall (default: pearson)
}
  \item{min.rows}{
  Minimum rows in "data" frame to proceed (default: 20)
}
  \item{Iterations}{
  Reboot Iterations (default: 1000)
}
  \item{subset1}{
  Subset of cols from cavity1 to iterate on (default: c(0)== ALL)
}
  \item{subset2}{
  Subset of cols from cavity2 to iterate on (default: c(0)== ALL)
}
  \item{errthresh}{
	Threshold error if there is enough data to calculate cor an pval for certain i and k (default: 0.0001)
}
  \item{distributions}{
	Output Distribution file  (default:NA)
}
}
 
\value{

	Returns a dataframe with five columns:
	bug1, bug2, cor, p.value, q.value
}
\author{
  Emma Schwager <emma.schwager@gmail.com>
}
\examples{
	mydata<-matrix(rlnorm(21*15,1,1),nrow=21)
	Calculation_Results <- ccrepe (frmOne=mydata,Iterations=5)
}
