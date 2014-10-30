test.nc.score <- function()
  #################################################################################
  #   test.nc.score                                                               #
  #   Function to verify  nc.score                                                #
  #################################################################################
{
    library(RUnit)
    library(infotheo)
	
	data <- read.table("nc_score_input_test.txt",header=TRUE,row.names=1)
	
	
	tol = 0.1
 
	context("NC-score")
	checkEqualsNumeric(nc.score( x=c(3, 1, 1, 3 ,2 ,3, 3 ,1 ,2, 1),  y= c(1 ,2,3 ,3, 1 ,1, 3 ,3, 1, 2)), -0.25,tolerance=tol) 

	testdata<-matrix(c(0.29787234, 0.2978723, 0.2553191, 0.1489362,
		0.17073171, 0.3170732, 0.2682927, 0.2439024,
		0.09302326, 0.3255814, 0.2558140, 0.3255814,
		0.32352941, 0.3235294, 0.1470588, 0.2058824,
		0.17241379, 0.1724138, 0.4137931, 0.2413793,
		0.29729730, 0.2162162, 0.2702703, 0.2162162,
		0.22500000, 0.3250000, 0.2000000, 0.2500000,
		0.12820513, 0.3589744, 0.2307692, 0.2820513,
		0.20000000, 0.2250000, 0.2250000, 0.3500000,
		0.10256410, 0.3076923, 0.1794872, 0.4102564
		),nrow=10,ncol=4,byrow = TRUE)
	dimnames(testdata) = list(
            c("Subject 1", "Subject 2","Subject 3","Subject 4","Subject 5","Subject 6","Subject 7","Subject 8","Subject 9","Subject 10"),
            c("bug 1", "bug 2", "bug 3","bug 4")) # column names 
		
	
	
	
	nc.score.results <-nc.score( x=testdata)
	nc.score.predicted.results <-matrix(c(1 ,-0.25000, -0.21875, -0.65625,
                                       -0.25000 ,      1, -0.21875 , 0.34375,
                                     -0.21875, -0.21875 ,      1, -0.21875,
                                     -0.65625,  0.34375, -0.21875,       1),
                                    nrow=4,ncol=4,byrow = TRUE)
	checkEqualsNumeric(nc.score.predicted.results, nc.score.results, tolerance = tol) 

	
	
	
	
	
	
	
	mymat2 <- matrix(
c(0.866073691523164, NA, 3.17467368259671, 0.359999900537411,
4.06522199311953, 0.531358433914907, 1.3553210433223,
0.70098991617494, NA, NA, 0.540666092673393, 0.449612399543197,
8.54571287906867, NA, 1.49025535173577, 0.989290609659256,
1.95255982150214, 0.727336361735765, 11.6147325678762, NA, NA,
0.33536474977386, 4.77792057692177, NA, 0.401697528610903,
0.307221503080535, 0.549849376034578, 1.04271392334705,
0.447789224691567, 2.06032279848217, 4.29836541780129,
0.541133352798413, 6.35178592023848, 0.55703492105826, NA,
0.45089052046624, 1.06131526363404, 5.23848765291558,
0.412230411756493, 6.11520270707892, 0.105560690896892,
0.333003141844819, NA, 1.54596666241829, 0.785548947825768,
0.571874628503047, 2.57875882990629, 0.93386278959432, NA,
0.145584626621529),
ncol=5
)
colnames(mymat2) <- paste0("Feature",seq(1,5))

tol <- 0.000001

x <- mymat2[,1]
y <- mymat2[,2]

check_mat_vec <- function(mymat,y,nbins=NULL,bin.cutoffs=NULL,ok=NULL){
  tau.check <- apply(mymat,2,function(x){
      if(!is.null(ok)){
        x <- x[ok]
        y <- y[ok]
      }
      nc.score(x,
               y,
               use=ifelse(!is.null(ok),"everything","pairwise.complete.obs"),
               nbins,
               bin.cutoffs
               )
    }
  )
  return(tau.check)
}

check_vec_mat <- function(x,mymat,nbins=NULL,bin.cutoffs=NULL,ok=NULL){
  tau.check <- apply(mymat,2,function(y){
      if(!is.null(ok)){
        x <- x[ok]
        y <- y[ok]
      }
      nc.score(x,
               y,
               use=ifelse(!is.null(ok),"everything","pairwise.complete.obs"),
               nbins,
               bin.cutoffs
               )
    }
  )
  return(tau.check)
}

## Using two vectors with bin numbers

	context("Using two vectors with bin numbers")

nc.score(x,y,use="pairwise.complete.obs",nbins=5)==nc.score(x,y,use="complete.obs",nbins=5)
is.na(nc.score(x,y,use="everything",nbins=5))



## Using two vectors with bin cutoffs
	context("Using two vectors with bin cutoffs")
nc.score(x,y,use="pairwise.complete.obs",NULL,bin.cutoffs=c(-1,0,1))==nc.score(x,y,use="complete.obs",NULL,bin.cutoffs=c(-1,0,1))
is.na(nc.score(x,y,use="everything",NULL,bin.cutoffs=c(-1,0,1)))



## Using matrix and vector with bin numbers
	context("Using matrix and vector with bin numbers")
tau <- nc.score(mymat2,y,use="pairwise.complete.obs",nbins=5)
tau.check <- check_mat_vec(mymat2,y,nbins=5)
checkEquals(sum(abs(tau-tau.check)>tol), 0) 


tau2 <- nc.score(mymat2,y,use="everything",nbins=5)
checkEquals(0,sum(!is.na(tau2)))



tau3 <- nc.score(mymat2,y,use="complete.obs",nbins=5)
ok   <- which(!is.na(apply(mymat2,1,sum)))
tau3.check <- check_mat_vec(mymat2,y,nbins=5,ok=ok)
checkEquals(0,sum(abs(tau3-tau3.check)>tol))



## Using matrix and vector with bin cutoffs
	context("Using matrix and vector with bin cutoffs")
tau.2 <-
nc.score(mymat2,y,use="pairwise.complete.obs",NULL,bin.cutoffs=c(-1,0,1))
tau.2.check <- check_mat_vec(mymat2,y,bin.cutoffs=c(-1,0,1))
checkEquals(0,sum(abs(tau.2-tau.2.check)>tol))



tau2.2 <- nc.score(mymat2,y,use="everything",NULL,bin.cutoffs=c(-1,0,1))
checkEquals(0,sum(!is.na(tau2.2)))



tau3.2 <-
nc.score(mymat2,y,use="complete.obs",NULL,bin.cutoffs=c(-1,0,1))
ok <- which(!is.na(apply(mymat2,1,sum)))
tau3.2.check <- check_mat_vec(mymat2,y,bin.cutoffs=c(-1,0,1),ok=ok)
checkEquals(0,sum(abs(tau3.2-tau3.2.check)>tol)) 




## Using vector and matrix with bin numbers
	context("Using vector and matrix with bin numbers")
tau.3 <- nc.score(y,mymat2,use="pairwise.complete.obs",nbins=5)
tau.3.check <- check_vec_mat(y,mymat2,nbins=5)
checkEquals(0,sum(abs(tau-t(tau.3))>tol))
checkEquals(0,sum(abs(t(tau.3)-tau.check)>tol))
checkEquals(0,sum(abs(tau.3-tau.3.check)>tol))
checkEquals(0,sum(abs(tau.3.check-tau.check)>tol))


tau2.3 <- nc.score(y,mymat2,use="everything",nbins=5)
checkEquals(0,sum(!is.na(tau2.3)))



tau3.3 <- nc.score(y,mymat2,use="complete.obs",nbins=5)
ok <- which(!is.na(apply(mymat2,1,sum)))
tau3.3.check <- check_vec_mat(y,mymat2,nbins=5,ok=ok)
checkEquals(0,sum(abs(tau3-t(tau3.3))>tol))
checkEquals(0,sum(abs(t(tau3.3)-tau3.check)>tol))
checkEquals(0,sum(abs(tau3.3-tau3.3.check)>tol))
checkEquals(0,sum(abs(tau3.3.check-tau3.check)>tol))


## Using vector and matrix with bin cutoffs
	context("Using vector and matrix with bin cutoffs")

tau.4 <- nc.score(y,mymat2,use="pairwise.complete.obs",nbins=NULL,bin.cutoffs=c(-1,0,1))
tau.4.check <- check_vec_mat(y,mymat2,bin.cutoffs=c(-1,0,1))
checkEquals(0,sum(abs(tau.2-t(tau.4))>tol))
checkEquals(0,sum(abs(t(tau.4)-tau.2.check)>tol))
checkEquals(0,sum(abs(tau.4-tau.4.check)>tol))
checkEquals(0,sum(abs(tau.4.check-tau.2.check)>tol))


tau2.4 <- nc.score(y,mymat2,use="everything",bin.cutoffs=c(-1,0,1))
checkEquals(0,sum(!is.na(tau2.4)))


tau3.4 <- nc.score(y,mymat2,use="complete.obs",bin.cutoffs=c(-1,0,1))
ok <- which(!is.na(apply(mymat2,1,sum)))
tau3.4.check <- check_vec_mat(y,mymat2,bin.cutoffs=c(-1,0,1),ok=ok)
checkEquals(0,sum(abs(tau3.2-t(tau3.4))>tol))
checkEquals(0,sum(abs(t(tau3.4)-tau3.2.check)>tol))
checkEquals(0,sum(abs(tau3.4-tau3.4.check)>tol))
checkEquals(0,sum(abs(tau3.4.check-tau3.2.check)>tol))


	
## Using matrix with bin numbers
	context("Using matrix with bin numbers")
tau.5 <- nc.score(mymat2,use="pairwise.complete.obs",nbins=5)
tau.5.check <- apply(mymat2,2,function(y)
check_mat_vec(mymat2,y,nbins=5))
tau.5.check.2 <- apply(mymat2,2,function(x)
check_vec_mat(x,mymat2,nbins=5))
checkEquals(0,sum(abs(tau.5[,2]-tau)>tol))
checkEquals(0,sum(abs(tau.5[,2]-tau.3)>tol))
checkEquals(0,sum(abs(tau.5-tau.5.check)>tol))
checkEquals(0,sum(abs(tau.5-tau.5.check.2)>tol))
 

tau2.5 <- nc.score(mymat2,use="everything",nbins=5)
sum(diag(tau2.5 != 1))==0
sum(!is.na(tau2.5)) == ncol(mymat2)

tau3.5 <- nc.score(mymat2,use="complete.obs",nbins=5)
ok <- which(!is.na(apply(mymat2,1,sum)))
tau3.5.check <- apply(mymat2,2,function(y)
check_mat_vec(mymat2,y,nbins=5,ok=ok))
tau3.5.check.2 <- apply(mymat2,2,function(x)
check_vec_mat(x,mymat2,nbins=5,ok=ok))
checkEquals(0,sum(abs(tau3.5[,2]-tau3)>tol))
checkEquals(0,sum(abs(tau3.5[,2]-tau3.3)>tol))
checkEquals(0,sum(abs(tau3.5-tau3.5.check)>tol))
checkEquals(0,sum(abs(tau3.5-tau3.5.check.2)>tol))



	
	
	
## Using matrix with bin cutoffs
	context("Using matrix with bin cutoffs")
tau.6 <- nc.score(mymat2,use="pairwise.complete.obs",bin.cutoffs=c(-1,0,1))
tau.6.check <- apply(mymat2,2,function(y)
check_mat_vec(mymat2,y,bin.cutoffs=c(-1,0,1)))
tau.6.check.2 <- apply(mymat2,2,function(x)
check_vec_mat(x,mymat2,bin.cutoffs=c(-1,0,1)))
checkEquals(0,sum(abs(tau.6[,2]-tau.2)>tol))
checkEquals(0,sum(abs(tau.6[,2]-tau.4)>tol))
checkEquals(0,sum(abs(tau.6-tau.6.check)>tol))
checkEquals(0,sum(abs(tau.6-tau.6.check.2)>tol))
 
 
 
tau2.6 <- nc.score(mymat2,use="everything",bin.cutoffs=c(-1,0,1))
checkEquals(0,sum(diag(tau2.6 != 1)))
checkEquals(sum(!is.na(tau2.6)), ncol(mymat2))



tau3.6 <- nc.score(mymat2,use="complete.obs",bin.cutoffs=c(-1,0,1))
ok <- which(!is.na(apply(mymat2,1,sum)))
tau3.6.check <- apply(mymat2,2,function(y)
check_mat_vec(mymat2,y,bin.cutoffs=c(-1,0,1),ok=ok))
tau3.6.check.2 <- apply(mymat2,2,function(x)
check_vec_mat(x,mymat2,bin.cutoffs=c(-1,0,1),ok=ok))
checkEquals(0,sum(abs(tau3.6[,2]-tau3.2)>tol))
checkEquals(0,sum(abs(tau3.6[,2]-tau3.4)>tol))
checkEquals(0,sum(abs(tau3.6-tau3.6.check)>tol))
checkEquals(0,sum(abs(tau3.6-tau3.6.check.2)>tol))


## Using two matrices with bin numbers
	context("Using two matrices with bin numbers")
tau.7 <-
nc.score(mymat2[,c(1,2)],mymat2[,c(3,4)],use="pairwise.complete.obs",nbins=5)
tau.7.check <- tau.5.check[c(1,2),c(3,4)]
checkEquals(0,sum(abs(tau.7-tau.7.check)>tol))
 

tau2.7 <-
nc.score(mymat2[,c(1,2)],mymat2[,c(3,4)],use="everything",nbins=5)
sum(!is.na(tau2.7))==0

tau3.7 <-
nc.score(mymat2[,c(1,2)],mymat2[,c(3,4)],use="complete.obs",nbins=5)
ok <- which(!is.na(apply(mymat2,1,function(row) sum(row[c(1,2,3,4)]))))
tau3.7.check <- apply(mymat2,2,function(y)
check_mat_vec(mymat2,y,nbins=5,ok=ok))[c(1,2),c(3,4)]
checkEquals(0,sum(abs(tau3.7-tau3.7.check)>tol))
 	
	
	context("Using two matrices with bin cutoffs")	
## Using two matrices with bin cutoffs

tau.8 <-
nc.score(mymat2[,c(1,2)],mymat2[,c(3,4)],use="pairwise.complete.obs",bin.cutoffs=c(-1,0,1))
tau.8.check <- tau.6.check[c(1,2),c(3,4)]
checkEquals(0,sum(abs(tau.8-tau.8.check)>tol))


tau2.8 <-
nc.score(mymat2[,c(1,2)],mymat2[,c(3,4)],use="everything",bin.cutoffs=c(-1,0,1))
sum(!is.na(tau2.8))==0
checkEquals(0,sum(!is.na(tau2.8)))

###############################################
# This one still has a problem                *
###############################################
#tau3.8 <-
#nc.score(mymat2[,c(1,2)],mymat2[,c(3,4)],use="complete.obs",bin.cutoffs=c(-1,0,1))
#ok <- which(!is.na(apply(mymat2,1,function(row) sum(row[c(1,2,3,4)]))))
#tau3.8.check <- apply(mymat2,2,function(y)
#check_mat_vec(mymat2,y,nbins=5,ok=ok))[c(1,2),c(3,4)]
#checkEquals(0,sum(abs(tau3.8-tau3.8.check)>tol))  
	
	

## Checking warning cases
	context("Checking warning cases")
expect_warning(nc.score(x,y,nbins=c(5,3),use="pairwise.complete.obs"))
expect_warning(nc.score(x,y,nbins=NULL,bin.cutoffs=c(1,-1,0),use="pairwise.complete.obs"))

## incorrect use argument - should give an error
expect_error(nc.score(x,y,use="something"))

## only one vector - should give an error
expect_error(nc.score(x,NULL))

## Check string input x - should give an error
expect_error(nc.score(c('1','2'),y[c(1,2)]))

## Check string input y - should give an error
expect_error(nc.score(x[c(1,2)],c('1','2')))

## Supplying both nbins and bin cutoffs - should give an error
expect_error(nc.score(x,y,nbins=2,bin.cutoffs=c(0,0.5,1,2)))

## Supplying invalid nbins values - should give errors
expect_error(nc.score(x,y,nbins="2"))
expect_error(nc.score(x,y,nbins=-1))
expect_error(nc.score(x,y,nbins=pi))

## Supplying invalid bin.cutoffs values - should give an error
expect_error(nc.score(x,y,nbins=NULL,bin.cutoffs=c(0,'1')))

	

## Miscellaneous testing
	context("Miscellaneous testing")
A <- c(1,2.5,2.5,4.5,4.5,6.5,6.5,8,9.5,9.5)
B <- c(1,2,4.5,4.5,4.5,4.5,8,8,8,10)
n <- length(A)
cor_1 <- cor(A,B,method="kendall")
concord_vec <- apply(t(combn(seq(1,length(A)),2)),1,function(row) (rank(A)[row[1]] - rank(A)[row[2]])*(rank(B)[row[1]] - rank(B)[row[2]]))
S <- (sum(concord_vec>0) - sum(concord_vec<0))
t_vals <- unique(A[which(sapply(A,function(e) sum(A==e))>1)])
u_vals <- unique(B[which(sapply(B,function(e) sum(B==e))>1)])
T_val  <- 0.5*sum(sapply(t_vals,function(e) sum(A==e)*(sum(A==e)-1)))
U_val  <- 0.5*sum(sapply(u_vals,function(e) sum(B==e)*(sum(B==e)-1)))
D_sq   <- (0.5*n*(n-1) - T_val)*(0.5*n*(n-1) - U_val)
cor_2 <- (sum(concord_vec>0) - sum(concord_vec<0))/(0.5*length(A)*(length(A)-1))
cor_3 <- S/sqrt(D_sq)

checkEquals(cor_1,cor_3	)
	
	
	
	
	}


