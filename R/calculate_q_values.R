calculate_q_values <-
function(CA, Gamma)
#**********************************************************************
#	Calculate Q values    				                              *
#**********************************************************************
{
	m = length(CA$p.values)					#m is the total number of p-values
	ln_m = log(m)									#log m
	ln_m_Plus_Gamma = ln_m + Gamma					
	SortedVector = sort(CA$p.values,index.return = TRUE)	#Sort the entire PValues matrix into a vector
	KVector = seq(1,m)						#A vector with the total number of entries in the PValues matrix
	QValues = SortedVector$x*m*ln_m_Plus_Gamma/KVector		#Calculate a vector containing the Q values
	QValuesArranged = rep(-1,m)
	QValuesArranged[SortedVector$ix] = QValues
	CA$q.values<-matrix(QValuesArranged, nrow= nrow(CA$p.values), byrow=TRUE)
	return(CA)
}
