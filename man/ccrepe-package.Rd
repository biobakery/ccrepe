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
 Calculate species-level co-variation and co-exclusion patterns based on an extension of the checkerboard score to ordinal data
}
\details{
\tabular{ll}{
Package: \tab ccrepe\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-04-18\cr
License: \tab MIT\cr
}
   An example of how to invoke ccrepe:\cr 
   mydata<-matrix(rlnorm(21*15,1,1),nrow=21)\cr
   ccrepe.results <-ccrepe(data1=mydata,iterations = 100)\cr
   An example of how to invoke nc.score:\cr
   nc.score.result <-nc.score (data1=mydata)\cr
}
\author{
Emma Schwager <emma.schwager@gmail.com>, Craig Bielski<craig.bielski@gmail.com>

Maintainer: Emma Schwager <emma.schwager@gmail.com>, Craig Bielski<craig.bielski@gmail.com>, George Weingart<george.weingart@gmail.com>
}
