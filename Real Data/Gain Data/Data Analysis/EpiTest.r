EpiTest = function(Y,X,geno1,geno2){

  N     	= nrow(Y)								# the number of individual

	Simi1	= TypicalIBS(geno1)
	S1		= Simi1$S
  H1    = Simi1$H
  R1    = Simi1$R
	q1		= Simi1$L

	Simi2	= TypicalIBS(geno2)
	S2		= Simi2$S
  H2    = Simi2$H
  R2    = Simi2$R
	q2		= Simi2$L
  
  #--------------------------------------------------------------------------
  # There are some concerns whether we need to use the trick to facilatating|
  # the inverse of V0: 1). When the sample size is to large, say, more than |
  # 1,000, just like the case in gain data, directly inverse V0 is OK but re|
  # -arly slow. In this case, we should prefer using the trick. 2). But usin|
  # -g the trick use a fatal problem is that S-S* where S*=HRH' may be far a|
  # -way from 0, So there I set a threshold, if the range of (S-S*) is blew |
  # e-10, then I think it is OK to use the trick, otherwise, I need to use t|
  # -he direct inverse anyway.                                              |
  #--------------------------------------------------------------------------
  
  FakeS1 = H1%*%R1%*%t(H1)
  diff1  = max(abs(range(FakeS1-S1)))
  
  FakeS2 = H2%*%R2%*%t(H2)
  diff2  = max(abs(range(FakeS2-S2)))
  
  TrickOK= (diff1<=1e-10) && (diff2<=1e-10) 
  
  if (TrickOK){
    cat("We can use the trick for this Genes combination","\n")
  }
  else{
    cat("diff1=",diff1,"diff2=",diff2,"so we can not use the trick","\n")
  }

	S12 	 = S1 * S2									# The similarity matrix for the interaction

	p 		 = ncol(X) 								# The number of fix effect

	P_X 	 = X %*% solve(crossprod(X,X))%*% t(X) 				# The projection of X

	#----------------------- Initialize Tau and sigma ----------------------------

	Tau_A_new = 1                        
	Tau_B_new = 1
	sigma_new = 1

	#------------------------------- updata --------------------------------------

	repeat{		
		Tau_A		  = Tau_A_new
		Tau_B 	  = Tau_B_new
		sigma 	  = sigma_new
		
		
    if (TrickOK) {
      R 			  = rbind(cbind(Tau_A*R1,array(0,dim <- c(q1,q2))),cbind(array(0,dim <- c(q2,q1)),Tau_B*R2))
		  H 			  = cbind(H1,H2)
		  inv_V0	  = 1/sigma*(diag(N)-1/sigma*H %*% R %*% solve(diag(q1+q2)+1/sigma*t(H)%*%H%*%R)%*%t(H))
    }
    else{
      V0        = Tau_A*S1 + Tau_B*S2 + sigma*diag(N)
      inv_V0    = solve(V0)  
    }
    
    inv_V0X 	= inv_V0%*%X                   #    To speed up the computation
		P 		  	= inv_V0-inv_V0X %*% solve(t(X)%*%inv_V0X)%*%t(inv_V0X)

    Tau_A_new = as.numeric(1/q1 * (Tau_A^2 * t(Y) %*% P %*% S1 %*% P %*% Y +  Tau_A*q1 - sum(diag(Tau_A^2 * P %*% S1))))
    Tau_B_new = as.numeric(1/q2 * (Tau_B^2 * t(Y) %*% P %*% S2 %*% P %*% Y +  Tau_B*q2 - sum(diag(Tau_B^2 * P %*% S2))))
    
    Y_star    = Y - Tau_A_new * S1 %*% P %*% Y - Tau_B_new * S2 %*% P %*% Y
    
    IP_X      = diag(N)-P_X
    VA        =Tau_A*S1+Tau_B*S2-(Tau_A*S1+Tau_B*S2) %*% P %*% (Tau_A*S1+Tau_B*S2)
    sigma_new = as.numeric(1/(N-p) * (t(Y_star) %*% IP_X %*% Y_star + sum(diag(IP_X%*%VA))))
    
#     cat(Tau_A_new,Tau_B_new,sigma_new,"\n")
    
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