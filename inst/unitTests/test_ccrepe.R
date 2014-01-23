test.ccrepe <- function()
{
	testdata<-matrix(c(
		0.29787234, 0.2978723, 0.2553191, 0.1489362,
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
		
		
	#***********************************************************
	#* ccrepe tests against testdata                           *
	#***********************************************************
	ccrepe.results <-  ccrepe  (x=testdata,   min.subj=10)
	tol = 0.15 
	p.values.results <- matrix(c(NA ,0.87723131, 0.5625562, 0.09091441,
							0.87723131, NA, 0.4573508, 0.07103074,
							0.56255622, 0.45735082,        NA, 0.90373585,
							0.09091441, 0.07103074, 0.9037359,         NA),
								nrow=4,ncol=4,byrow = TRUE)
	
	q.values.results <- matrix(c(NA, 2.472189, 2.066321, 0.8325509,
								2.4721889, NA, 2.273340, 1.1677655,	
								2.0663210, 2.273340,  NA, 2.2336444,
								0.8325509, 1.167766, 2.233644,  NA), 
								nrow=4,ncol=4,byrow = TRUE)
								
	sim.score.results <- matrix(c(1.0000000, -0.3454545, -0.1636364, -0.7818182,
								-0.3454545,  1.0000000, -0.4424242,  0.2727273,
								-0.1636364, -0.4424242 , 1.0000000, -0.2727273,
								-0.7818182,  0.2727273, -0.2727273,  1.0000000),
								nrow=4,ncol=4,byrow = TRUE)
				
	z.stat.results <-matrix(c(   NA,  0.1544800,  0.5790488, -1.6905938,
								0.1544800,   NA, -0.7432162,  1.8052809,
								0.5790488, -0.7432162,         NA,  0.1209434,
								-1.6905938,  1.8052809,  0.1209434,         NA),
								nrow=4,ncol=4,byrow = TRUE)							


	checkEqualsNumeric(p.values.results,ccrepe.results$p.values, tolerance = tol)
	checkEqualsNumeric(q.values.results,ccrepe.results$q.values, tolerance = tol)
	checkEqualsNumeric(sim.score.results,ccrepe.results$sim.score, tolerance = tol)
	checkEqualsNumeric(z.stat.results,ccrepe.results$z.stat, tolerance = tol)
	
	

	
	
	
	#***********************************************************
	#* nc.score tests                                          *
	#***********************************************************
	checkEqualsNumeric(-0.3809524, nc.score( x=c(3, 1, 1, 3 ,2 ,3, 3 ,1 ,2, 1),  y= c(1 ,2,3 ,3, 1 ,1, 3 ,3, 1, 2)), tolerance = 0.0001)

	nc.score.results <-nc.score( x=testdata)
	nc.score.predicted.results <-matrix(c(NA, -0.3809524, -0.3333333, -0.8354978,
			-0.3809524 , NA ,-0.3333333,  0.3593074,
			-0.3333333, -0.3333333, NA, -0.3333333,
			-0.8354978,  0.3593074, -0.3333333,  NA),
			nrow=4,ncol=4,byrow = TRUE)
 
	checkEqualsNumeric(nc.score.predicted.results,nc.score.results, tolerance = tol)
 
	
}

context("Permutation p-values")

test_that("Renormalization gives normalized data without missing", {
    data <- matrix(1:10,nrow=2,byrow=TRUE)
    data.norm <- ccrepe_norm(data)

    expect_that( is.matrix(data.norm), is_true())
    expect_that( apply(data.norm,1,sum), equals(c(1,1)) )
    expect_that( data.norm, equals( matrix(c(1/15,2/15,3/15,4/15,5/15,6/40,7/40,8/40,9/40,0.25),
                                           byrow=TRUE,
                                           nrow=2) ) )
})

test_that("Renormalization works with missing data", {
    data <- matrix(1:10,nrow=2,byrow=TRUE)
    data[1,2] = NA
    data.norm <- ccrepe_norm(data)

    expect_that( data.norm[1,], equals(rep(0,ncol(data))) )
})


test_that("Permutation of one matrix permutes the data by columns", {
    data <- matrix(1:12, ncol=3)
    permute.id.matrix <- matrix(c(4,3,2,1,2,1,4,3,2,3,1,4),ncol=3)
    data.permute <- permute(data,permute.id.matrix)

    expect_that( is.matrix(data.permute), is_true() )
    expect_that( data.permute, equals(matrix(c(4,3,2,1,6,5,8,7,10,11,9,12),ncol=3)) )
})

v_dist.na   <- c(3.4,2,NA,2,4.6,0,10,2,8,NA,-2,0,NA)
v_dist      <- c(3.4,2,2,4.6,0,10,2,8,-2,0)
v_dist.null <- NA
obs.value1 <- 3.4             # high p-value, in dist
obs.value2 <- 10              # low p-value, in dist
obs.value3 <- 4               # high p-value, not in dist
obs.value4 <- -5              # low p-value, not in dist
dist.value1 <- 2              # repeated value
dist.value2 <- 0              # repeated value

test_that("get_count returns a correct value with no missing", {

    count1 <- get_count(dist.value1, v_dist)
    count2 <- get_count(dist.value2, v_dist)
    count3 <- get_count(obs.value3, v_dist)
    count4 <- get_count(obs.value1, v_dist)

    expect_that( count1, equals(3) )
    expect_that( count2, equals(2) )
    expect_that( count3, equals(0) )
    expect_that( count4, equals(1) )
})

test_that("get_perm_p.value works with no missing", {

    p.value1 <- get_perm_p.value(v_dist,obs.value1)
    p.value2 <- get_perm_p.value(v_dist,obs.value2)
    p.value3 <- get_perm_p.value(v_dist,obs.value3)
    p.value4 <- get_perm_p.value(v_dist,obs.value4)

    expect_that( p.value1, equals(1) )
    expect_that( p.value2, equals(.1) )
    expect_that( p.value3, equals(.9) )
    expect_that( p.value4, equals(0) )
})

test_that("get_perm_p.value works with missing", {

    p.value1 <- get_perm_p.value(v_dist.na,obs.value1)
    p.value2 <- get_perm_p.value(v_dist.na,obs.value2)
    p.value3 <- get_perm_p.value(v_dist.na,obs.value3)
    p.value4 <- get_perm_p.value(v_dist.na,obs.value4)

    expect_that( p.value1, equals(1) )
    expect_that( p.value2, equals(.1) )
    expect_that( p.value3, equals(.9) )
    expect_that( p.value4, equals(0) )
})

test_that("get_perm_p.value returns missing if only missing in v_dist", {

    p.value <- get_perm_p.value(v_dist.null,obs.value1)

    expect_that( is.na(p.value), is_true() )
})

test_that("get_perm_p.value returns missing if obs.value is missing", {

    p.value <- get_perm_p.value(v_dist, NA)

    expect_that( is.na(p.value), is_true() )
})

 
