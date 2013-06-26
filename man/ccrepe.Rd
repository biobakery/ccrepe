\name{ccrepe}
\alias{ccrepe}
\title{
  ccrepe 
   When method=cor: A function to calculate correlation, p-values and q-values of the relative abundance of bugs in one or two cavities.
   When method=nc.score: A function to calculate species-level co-variation and co-exclusion patterns based on an extension of the checkerboard score to ordinal data  
}
\description{
   ccrepe uploads the relative abundances of bugs in one or two cavities and calculates their correlation, p-value and q-value using permutation\cr
   and rebooting of the data\cr
 }
\usage{
 ccrepe(
 x = NA,
 y = NA,
 method = cor,
 method.args = list(),
 min.subj = 20,
 iterations = 1000,
 errthresh = 1e-04,
 verbose = FALSE,
 iterations.gap = 100,
 distributions = NA)
 
}
\arguments{
  \item{x}{
   First dataframe containing the relative abundances in cavity1 : Columns are bugs, rows are samples\cr
   The first column is considered data\cr 
}
  \item{y}{
   Second dataframe (optional) containing the relative abundances in cavity2: Columns are bugs, rows are samples \cr
   The first column is considered data 
   This dataframe will be merged with the first dataframe by Row Number, therefore same rows must pertain to same subjects
   The matrices are merged by row name (Which defaults to row number)  
}
 
  \item{method}{
   Measurement function, such as cor\cr
   some properties of the function:\cr
   1. Be able to take either two inputs which are vectors or one input which is either a matrix or a data frame\cr
   2.In the case of two inputs, return a single number\cr
   3.In the case of one input, return a matrix in which the (i,j)th entry is the similarity score for column i and column j in the original matrix\cr
   4.Resulting matrix must be symmetric\cr
   5.The inputs must be named x and y\cr
}
  \item{method.args}{
	A list of arguments for the measurement function.
	For example:  In the case of cor, the following would be acceptable:
	method.args = list(method='spearman',use='complete.obs'	)	
}
  \item{min.subj}{
  Minimum number of subjects that must be present (default: 20)
}
  \item{iterations}{
  Reboot iterations (default: 1000)
}
  \item{errthresh}{
	errThresh   (default= 0.0001)  Represents the probability of getting all 0's in a given bootstrapped column.
	If a bug/feature/column has too many 0's, its probability of sampling all 0's will be too large; it is therefore not used
}
  \item{verbose}{
	Indicator whether the User requested verbose Output  (default:FALSE)
}
 \item{iterations.gap}{
	If output is verbose - number of iterations after issue a status message (Default=100 - displayed only if verbose=TRUE)
}
  \item{distributions}{
	Output Distribution file  (default:NA)
}
}
 
\value{
    For ccrepe:\cr
	Returns a dataframe containing the calculations results and the parameters used.\cr
	Parameters shown:\cr
		min.subj, iterations, errthresh, method, subset.cols.1  (and subset.cols.2 -if applicable)\cr
		iterations.gap and method.args.method('spearman','kendall' or 'pearson')\cr
		Their descriptions are above\cr
	The calculation results consist of:\cr
		cor.matrix, p-values.matrix,  q-values.matrix  
	\cr
	\cr
	For nc.score\cr
	It returns a dataframe containing the calculation results and the parameters used.\cr
	Parameters shown:\cr
	min.abundance, min.samples and adBons\cr
	In addition, the function displays the matrices of the filtered data x.filtered and y.filtered\cr
	Filtering consists of eliminating NA's as well as observations that do not meet \cr
	the minimum abundance threshold: Selection of species with abundance has to be  >= min.abundance \cr
	in more than min.samples percent of samples\cr
	The calculation results consists of matrix or vectors of normalized scores\cr
		
}
\author{
  Emma Schwager <emma.schwager@gmail.com> 
}
\examples{
	data1<-matrix(rlnorm(21*23,1,1),nrow=21)
	data2<-matrix(rlnorm(23*24,1,1),nrow=24)
	ccrepe.results <-ccrepe  (x=data1,iterations=100) 
}
