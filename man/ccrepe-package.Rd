\name{ccrepe-package}
\alias{ccrepe-package}
\docType{package}
\title{
 The package has two functions 
 1. ccrepe: Calculate correlations, p-values and q-values for relative abundances of bugs in one or two body sites
 2. nc.score: Calculate Species-Level Co-Variation and Co-Exclusion Patterns in the Human Microbiome
}
\description{
 Function number 1 (ccrepe):    Calculate correlations, p-values and q-values for relative abundances of bugs in one or two body sites using bootstrap and permutation matrices of the data
 Function number 2 (nc.score):  Calculate species-level co-variation and co-exclusion patterns based on an extension of the checkerboard score to ordinal data
}
\details{
\tabular{ll}{
Package: \tab ccrepe\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-04-18\cr
License: \tab MIT\cr
}
   An example of how to invoke ccrepe:
   mdPvalues <- ccrepe (frmOne=data3,ProcessType='ccrepe')
   An example of how to invoke nc.score:
   Correlation.Matrix <- ccrepe(frmOne=data3,ProcessType='nc.score') 
}
\author{
Emma Schwager <emma.schwager@gmail.com>, Craig Bielski<craig.bielski@gmail.com>

Maintainer: Emma Schwager <emma.schwager@gmail.com>, Craig Bielski<craig.bielski@gmail.com>, George Weingart<george.weingart@gmail.com>
}
