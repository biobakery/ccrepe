\name{ccrepeSampleTestFunction}
\alias{ccrepeSampleTestFunction}
\title{
 ccrepeSampleTestFunction - Simple example of a test measurent function to be used with ccrepe       
}
\description{
This  simple example of a test measurent function to be used with ccrepe\cr       
 used in the same fashion that cor would be used \cr                                  
 Some properties of the function:  \cr                                                
 1. Be able to take either two inputs which are vectors or one input which  \cr       
 is either a matrix or a data frame  \cr                                            
 3.In the case of one input, return a matrix in which the (i,j)th entry  \cr          
   is the similarity score for column i and column j in the original matrix \cr        
 4.Resulting matrix must be symmetric \cr                                              
 5.The inputs must be named x and y \cr   
}
\usage{
ccrepeSampleTestFunction(x, y = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 x is a vector or a matrix
}
  \item{y}{
 y is a vector.\cr
 if y selected then x must be a vector too
}
}
\value{
 If x and y are vectors it returns a number: 0.5
 If x is a matrix it returns a matrix of all 0.5
}

\author{
  Emma Schwager <emma.schwager@gmail.com> 
}
