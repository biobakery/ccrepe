\name{ccrepe-package}
\alias{ccrepe-package}
\docType{package}
\title{
 The package has two functions\cr 
 1. ccrepe: Calculate correlations, p-values and q-values for relative abundances of bugs in one or two body sites \cr 
 2. nc.score: Calculate Species-Level Co-Variation and Co-Exclusion Patterns in the Human Microbiome
}
\description{
 Function number 1 (ccrepe):\cr 
     Calculate correlations, p-values and q-values for relative abundances of bugs in one or two body sites using bootstrap and permutation matrices of the data.\cr 
 Function number 2 (nc.score):\cr 
     Calculate species-level co-variation and co-exclusion patterns based on an extension of the checkerboard score to ordinal data\cr
 The calculations of ccrepe or nc.score are controlled via the parameter method\cr
 When method=cor the package performs ccrepe calculations
 When method=nc.score, the package perform nc.score calculations
}
\details{
\tabular{ll}{
Package: \tab ccrepe\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-04-18\cr
License: \tab MIT\cr
}
   An example of how to invoke the ccrepe calculations:\cr 
    data1<-matrix(rlnorm(21*5,1,1),nrow=21)\cr
	ccrepe.results <-ccrepe  (x=data1, method = cor ,iterations=100)\cr
	\cr
   An example of how to invoke the nc.score calculations:\cr 	
   \cr
    data1<-matrix(rlnorm(21*5,1,1),nrow=21)\cr
	data2<-matrix(rlnorm(5*17,1,1),nrow=17)\cr
	nc.score.results <- nc.score( x=data1, y=data2 )
}
\author{
Emma Schwager <emma.schwager@gmail.com>, Craig Bielski<craig.bielski@gmail.com>\cr
Maintainer: Emma Schwager <emma.schwager@gmail.com>, Craig Bielski<craig.bielski@gmail.com>, George Weingart<george.weingart@gmail.com>
}
