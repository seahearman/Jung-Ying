EpiTest = function(Y,X,geno1,geno2){

  N     	= nrow(Y)								# the number of individual

	Simi1	= TypicalIBS(geno1)
	S1		= Simi1$S
	q1		= Simi1$L

	Simi2	= TypicalIBS(geno2)
	S2		= Simi2$S
	q2		= Simi2$L

	S12 	= S1 * S2									# The similarity matrix for the interaction

	p 		= ncol(X) 								# The number of fix effect

	P_X 	= X %*% solve(crossprod(X,X))%*% t(X) 				# The projection of X

	#----------------------- Initialize Tau and sigma ----------------------------

	Tau_A_new = 1                        
	Tau_B_new = 1
	sigma_new = 1

	#------------------------------- updata --------------------------------------

	repeat{		
		Tau_A		  = Tau_A_new
		Tau_B 	  = Tau_B_new
		sigma 	  = sigma_new
		
		V0        = Tau_A*S1 + Tau_B*S2 + sigma*diag(N)
    inv_V0    = solve(V0)
    inv_V0X 	= inv_V0%*%X                   #    To speed up the computation
		P 		  	= inv_V0-inv_V0X %*% solve(t(X)%*%inv_V0X)%*%t(inv_V0X)

    Tau_A_new = as.numeric(1/q1 * (Tau_A^2 * t(Y) %*% P %*% S1 %*% P %*% Y +  Tau_A*q1 - sum(diag(Tau_A^2 * P %*% S1))))
    Tau_B_new = as.numeric(1/q2 * (Tau_B^2 * t(Y) %*% P %*% S2 %*% P %*% Y +  Tau_B*q2 - sum(diag(Tau_B^2 * P %*% S2))))
    
    Y_star    = Y - Tau_A_new * S1 %*% P %*% Y - Tau_B_new * S2 %*% P %*% Y
    
    IP_X      = diag(N)-P_X
    VA        =Tau_A*S1+Tau_B*S2-(Tau_A*S1+Tau_B*S2) %*% P %*% (Tau_A*S1+Tau_B*S2)
    sigma_new = as.numeric(1/(N-p) * (t(Y_star) %*% IP_X %*% Y_star + sum(diag(IP_X%*%VA))))
    
		diff_A = abs((Tau_A_new-Tau_A)/1)		# use the relative difference rather than the absoult difference
		diff_B = abs((Tau_B_new-Tau_B)/1)
		diff_s = abs((sigma_new-sigma)/1)
		if ((diff_A<0.00001) & (diff_B<0.00001) & (diff_s<0.00001)) break
	}

  P0 					  = P 
	T0 					  = 1/2*t(Y) %*% P0 %*% S12 %*% P0 %*% Y	#Get T0 under null hypothesis

  e						  = eigen(V0, symmetric=TRUE)
	V_eigen				= e$vectors
	V_square_root	= V_eigen %*% diag(sqrt(e$values)) %*% t(V_eigen)


	Weights_all		= eigen(1/2*V_square_root %*% P0 %*% S12 %*% P0 %*% V_square_root, symmetric=TRUE, only.value=TRUE)
	temp					= Weights_all$values
	temp2					= sort(temp,decreasing=TRUE)
	dim(temp2)		= c(N,1)
	big_enough 		= sum(temp>10^-3)			# Get the number of big eigen values. here, the threshold for "big" is 10^-3 
	Weights				= array(temp2[1:big_enough,1],dim=c(big_enough,1))

  p = liu(T0, Weights, h = rep(1, length(Weights)), delta = rep(0, length(Weights)))

	return(P=p)
	
}