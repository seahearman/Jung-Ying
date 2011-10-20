JointTest = function(Y,X,Simi1,Simi2){
  N 			= nrow(Y)															# the number of individual

	S1			= Simi1$S

	S2			= Simi2$S

	S12 		= S1 * S2														# The similarity matrix for the interaction
	Q 		= diag(N) - X %*% solve(t(X) %*% X) %*% t(X)

	sigma 	= as.numeric (t(Y) %*% Q %*% Y/(N-ncol(X)))

	P0 		= 1 / sigma * Q

	T0 		= 1 / (2*sigma^2) * t(Y) %*% Q %*% (S1 + S2 + S12) %*% Q %*% Y
	Weights_all = eigen(1 / (2 * sigma) * Q %*% (S1 + S2 + S12) %*% Q, symmetric=TRUE, only.value=TRUE)
	temp 	= Weights_all$values
	temp2 	= sort(temp,decreasing=TRUE)
	dim(temp2) = c(N,1)
	big_enough = sum(temp>10^-3)			# Get the number of big eigen values. here, the threshold for "big" is 10^-3 

	Weights = array(temp2[1:big_enough,1],dim=c(big_enough,1))

  p = liu(T0, Weights, h = rep(1, length(Weights)), delta = rep(0, length(Weights)))
	return(P=p)

}