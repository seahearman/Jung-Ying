OneGeneTest = function(Y,X,Simi1,Simi2){

  N   		= nrow(Y)															# the number of individual

	S1		= Simi1$S
	q1		= Simi1$L
  
	S2		= Simi2$S

	p 		= ncol(X) 														# The number of fix effect

	P_X 	= X %*% solve(crossprod(X,X))%*% t(X) 				# The projection of X

	#----------------------- Initialize Tau and sigma ----------------------------

	Tau_A_new = 0.5 
  sigma_new = 0.5

	#------------------------------- updata --------------------------------------

	repeat{		
		Tau_A		  = Tau_A_new
		sigma 	  = sigma_new
		
    V0        = Tau_A * S1 + sigma * diag(N)
    inv_V0    = solve(V0)
		inv_V0X 	= inv_V0%*%X                   #    To speed up the computation
		P 		  	= inv_V0-inv_V0X %*% solve(t(X)%*%inv_V0X)%*%t(inv_V0X)

    Tau_A_new = as.numeric(1/q1 * (Tau_A^2 * t(Y) %*% P %*% S1 %*% P %*% Y + Tau_A*q1 - sum(diag( Tau_A^2 * P %*% S1))))
    
    Y_star    = Y - Tau_A * S1 %*% P %*% Y
    
    IP_X      = diag(N)-P_X
    VA        =Tau_A*S1-(Tau_A*S1) %*% P %*% (Tau_A*S1)
    sigma_new = as.numeric(1/(N-p) * (t(Y_star) %*% IP_X %*% Y_star + sum(diag(IP_X%*%VA))))

		diff_A = abs((Tau_A_new-Tau_A)/Tau_A_new+1e-10)		# use the relative difference rather than the absoult difference
		diff_s = abs((sigma_new-sigma)/sigma_new+1e-10)
		if ((diff_A<0.0001)  & (diff_s<0.0001)) break 
    cat(diff_A,"\n")
 	}

	P0 					  = P 

	T0 					  = 1/2*t(Y) %*% P0 %*% S2 %*% P0 %*% Y	#Get T0 under null hypothesis
	e						  = eigen(V0, symmetric=TRUE)
	V_eigen				= e$vectors
	V_square_root	= V_eigen %*% diag(sqrt(e$values)) %*% t(V_eigen)


	Weights_all		= eigen(1/2*V_square_root %*% P0 %*% S2 %*% P0 %*% V_square_root, symmetric=TRUE, only.value=TRUE)
	temp					= Weights_all$values
	temp2					= sort(temp,decreasing=TRUE)
	dim(temp2)		= c(N,1)
	big_enough 		= sum(temp>10^-3)			# Get the number of big eigen values. here, the threshold for "big" is 10^-3 
	Weights				= array(temp2[1:big_enough,1],dim=c(big_enough,1))

  p = liu(T0, Weights, h = rep(1, length(Weights)), delta = rep(0, length(Weights)))
	return(P=p)
	
}