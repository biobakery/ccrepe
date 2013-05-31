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
   nc.score calculates species-level co-variation and co-exclusion patterns of two inputs based on an extension of the checkerboard score to ordinal data
 }
\usage{
 ccrepe(
 x = NA,
 y = NA,
 method = cor,
 method.args = list(),
 min.subj = 20,
 iterations = 1000,
 subset.cols.1 = c(0),
 subset.cols.2 = c(0),
 errthresh = 1e-04,
 verbose = FALSE,
 iterations.gap = 100,
 distributions = NA,
 bins = NA,
 min.abundance = 1e-04,
 min.samples = 0.1)
 
}
\arguments{
  \item{x}{
   First dataframe containing the relative abundances in cavity1 : Columns are bugs, rows are samples
}
  \item{y}{
   Second dataframe containing the relative abundances in cavity2: Columns are bugs, rows are samples
}
 
  \item{method}{
	If method=cor   the function invokes ccrepe calculations
	If method=nc.score the function invokes nc.score calculations
	If it is another function, provided bu the User, it must have the following characteristics:
	Needs to be able to take two vectors or one matrix; Does not have to be able to take two matrices
}
  \item{method.args}{
	A list of arguments for the cor function.
	Used to select the correlation method 
	Could be
	list(method = "pearson")
	list(method = "spearman")
	list(method = "kendall")
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
  \item{verbose}{
	Indicator whether the User requested verbose Output  (default:FALSE)
}
 \item{iterations.gap}{
	If output is verbose - number of iterations after issue a status message (Default=100 - displayed only if verbose=TRUE)
}
  \item{distributions}{
	Output Distribution file  (default:NA)
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
    For ccrepe:\cr
	Returns a dataframe containing the calculations results and the parameters used.\cr
	Parameters shown:\cr
		min.subj, iterations, errthresh, method, subset.cols.1  (and subset.cols.2 -if applicable)\cr
		iterations.gap and method.args.method('spearman','kendall' or 'pearson')\cr
		Their descriptions are above\cr
	The calculation results consist of:\cr
		PValues.matrix,  Q.matrix and  data.cor that contains:  bug1 bug2 cor   p.value  q.value\cr
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
  Emma Schwager <emma.schwager@gmail.com>,Craig Bielski<craig.bielski@gmail.com>
}
\examples{
	data1<-matrix(rlnorm(21*5,1,1),nrow=21)
	data2<-matrix(rlnorm(5*17,1,1),nrow=17)
	ccrepe.results <-ccrepe  (x=data1, method = cor ,iterations=100) 
	nc.score.results <- ccrepe( x=data1, y=data2, method = nc.score )
}
