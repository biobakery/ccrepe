cor.CalcFunction <-
function(x,y,method){
#*************************************************************************************
#*	Correlation function                                                         *
#*************************************************************************************
   	if (!is.matrix(x))
		{
		return (0)
		}
	if (missing(y))
		{
		return(cor(x, use = "complete.obs", method = method))		#For matrices
		}
	else
		{
		return(cor(x, y, use = "complete.obs", method = method))	#Between two elements
		}
	}
