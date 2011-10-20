# Xin Wang, bioinformatics, NCSU, Aug, 2011 I use this to have a try

library(stats)
library(MASS)
library(CompQuadForm)
rm(list= ls())
set.seed(4)

#--------------------------------------------------------------------------------------------------
#														FUNCTIONS																		|
#--------------------------------------------------------------------------------------------------

#-------------------- Function : weighted chi-square distribution---------------------
# The main idea for for this function is:
# 1. generate 8*10,000 chi square distributed random variables: Matrix C(this should be
#    outside of the function to save time, you don't need to generate lots of random va
#    -irables each time)
# 2. Acoording to the eigenvalues, calculated 10,000 weighted chi-square distributed ra
#    -ndom variables
# 3. given the Ti statistics, get the p-value.

weighted_chi2 = function(eigenvalues,Ti){

	c1 	= sum(eigenvalues)
	c2 	= sum(eigenvalues^2)
	c3 	= sum(eigenvalues^3)
	h 	= c2^3/c3^2
	y 	= (Ti-c1)*(h/c2)^(0.5)+h
	p 	= pchisq(y,df=h)
	return(p_value=p)
	
	C 	= rchisq(n=l*5000,df=1)
	dim(C) 	= c(5000,l)	# so now the dimension of C is L*5000

	weightC = C %*% Weights 
	N_w 		 = nrow(weightC)
	p 			= sum(as.numeric(T)<=weightC)/N_w
	
}

#-------------- Function : generate the similarity matrix from the genotype-----------
# This funcion is used to calculate the typical IBS from the given genotype

TypicalIBS = function (geno){
  #- if the genotype has more than 2 allelic data, we may want to calculated in other way
  #- Sepetating the genotype into 2 haplotype matrix ----

# 	hA 	  		= geno				# hA is used to how many A's in the allele, it equals to genotype
# 	ha 	  		= 2-geno			# ha is used to how many a's in the allele, it equals to 2 - genotype
# 	
# 	N_gene  = ncol(geno)    # the number of markers
# 	S_temp  = hA %*% t(hA) + ha %*% t(ha)
# 	S         	=  S_temp/(4*N_gene)
  
  #--- how to generate the S matrix when the genotype is multi-allelic
  
  N = nrow(geno)
  M = ncol(geno)
  
  S = array(0,dim=c(N,N))
  
  for (i in 1:N){
    for (j in 1:N){
      
      #- compute S[i,j]
      p1 = 1
      p2 = p1+1
      
      while (p2<=M){
        if ( (geno[i,p1]==geno[i,p2]) && (geno[i,p1]==geno[j,p1]) && (geno[j,p1]==geno[j,p2]))
          S[i,j] = S[i,j] + 2/M
        else if ((geno[i,p1]!=geno[j,p1]) && (geno[i,p1]!=geno[j,p2]) && (geno[i,p2]!=geno[j,p1]) && (geno[i,p2]!=geno[j,p2]))
          S[i,j] = S[i,j]
        else S[i,j] = S[i,j] + 1/M
        p1 = p2+1
        p2 = p1+1
      }
      if (i==j){
        S[i,j]=1
      }
    }
      
  }
	
  eg     =  eigen(S, symmetric=T)
  evalue =  eg$values
  le     =  evalue[1] # Biggest eigenvalue.
	if (le == 0){
		le 		= 1
	}
	tmpkey = (evalue / le) > 1e-7
	ll  	 		= sum(tmpkey) # Rank of SSS
	RRRR	 	= diag(evalue[tmpkey])
	HHHH   	= eg$vectors[,tmpkey]
	return(list(S = S, R = RRRR, H = HHHH, L = ll))
	
}

ModifiedIBS = function (geno){

	#- Sepetating the genotype into 2 haplotype matrix ----

	hA 	 		 = geno				# hA is used to how many A's in the allele, it equals to genotype
	ha 	  		= 2-geno			# ha is used to how many a's in the allele, it equals to 2 - genotype
	
	N 		  	= nrow(geno)    # sample size
	N_gene = ncol(geno)    # the number of markers
	
	S 		  	= array(0,dim=c(N,N))
	
	#--- Generate the Sij -----
	
	for (i in 1:N){
		for(j in i:N){
			
			tempS 				= (hA[i,]*hA[j,]+ha[i,]*ha[j,])/(4*N_gene)
			dim(tempS)  	=  c(1,N_gene)
			S[i,j] 				= sum(tempS)
			
			for (i1 in 1:(N_gene-1)){
				for (j1 in (i1+1):N_gene){
				
					S[i,j]		= S[i,j]+tempS[1,i1]*tempS[1,j1]
					
				}
			}
			
			S[j,i] 				= S[i,j]
		}
	}
    eg        =  eigen(S, symmetric=T)
    evalue  =  eg$values
    le        =  evalue[1] # Biggest eigenvalue.
	if (le == 0){
		le = 1
	}
	tmpkey = (evalue / le) > 1e-7
	ll  	 = sum(tmpkey) # Rank of SSS
	RRRR	 = diag(evalue[tmpkey])
	HHHH   = eg$vectors[,tmpkey]
	return(list(S = S, R = RRRR, H = HHHH, L = ll))
	
}


MaxIBS		 = function (geno){

	N_gene = ncol(geno)
	N 			= nrow(geno)
	
	h1A_ = array(0,dim=c(N,N_gene))			# h1A_ is used to indicate whether the 1st allele is A for gene
	h1_a = array(0,dim=c(N,N_gene))			# h1_a is used to indicate whether the 1st allele is a for gene
	h2A_ = array(0,dim=c(N,N_gene))			# h2A_ is used to indicate whether the 2nd allele is A for gene
	h2_a = array(0,dim=c(N,N_gene))			# h2_a is used to indicate whether the 2nd allele is a for gene

	h1_a = (geno>0)*1						# This is a small algorithm to generate the h1A_~h2_a from the genotype
	h2_a = geno-h1_a

	h1A_ = 1-h1_a
	h2A_ = 1-h2_a


	#----------------- Get the allele freq for each locus in each two gene-------------
	qA 	= array(1,dim=c(1,N_gene))		# calculate the freq of main  allele for gene
	qa	= array(1,dim=c(1,N_gene))		# calculate the freq of minor allele for gene



	wA 	=qA^(-0.75)			# using the freq^(-3/4) as the weight
	wa 	=qa^(-0.75)

	S_original = array(0,dim=c(N,N))						# S by {0, 0.5, 1} value set


	# Following is small loop to generate the similarity matrix for gene1 and gene2	

	for (i in 1:N_gene){
		temp1 = wA[i]*h1A_[,i] %*% t(h1A_[,i])+wa[i]*h1_a[,i] %*% t(h1_a[,i]) + wA[i]*h2A_[,i] %*% t(h2A_[,i])+wa[i]*h2_a[,i] %*% t(h2_a[,i])
		temp2 = wA[i]*h1A_[,i] %*% t(h2A_[,i])+wa[i]*h1_a[,i] %*% t(h2_a[,i]) + wA[i]*h2A_[,i] %*% t(h1A_[,i])+wa[i]*h2_a[,i] %*% t(h1_a[,i])
		S_original = S_original+ temp1*(temp1>=temp2)+temp2*(temp1<temp2)
	}

	S = S_original/(2*N_gene)
	
    eg        =  eigen(S, symmetric=T)
    evalue  =  eg$values
    le        =  evalue[1] # Biggest eigenvalue.
	if (le == 0){
		le = 1
	}
	tmpkey = (evalue / le) > 1e-7
	ll  	 = sum(tmpkey) # Rank of SSS
	RRRR	 = diag(evalue[tmpkey])
	HHHH   = eg$vectors[,tmpkey]
	return(list(S = S, R = RRRR, H = HHHH, L = ll))	

}

EpiTest = function(Y,X,geno1,geno2){


	return(P=p)
	
}
OneGeneTest = function(Y,X,geno1,geno2){

  N 			= nrow(Y)															# the number of individual

	Simi1	= TypicalIBS(geno1)
	S1		= Simi1$S
	H1		= Simi1$H 
	R1		= Simi1$R
	q1		= Simi1$L
  
  Simi2	= TypicalIBS(geno2)
	S2		= Simi2$S

	p 		= ncol(X) 														# The number of fix effect

	P_X 	= X %*% solve(crossprod(X,X))%*% t(X) 				# The projection of X

	#----------------------- Initialize Tau and sigma ----------------------------

	Tau_A_new = 1                        
	sigma_new = 1
#   S1 = H1 %*% R1 %*% t(H1)
	#------------------------------- updata --------------------------------------

	repeat{		
		Tau_A		  = Tau_A_new
		sigma 	  = sigma_new
		
		R 			  = Tau_A * R1
		H 			  = H1
		V0			  = H %*% R %*% t(H) + sigma * diag(N)

		inv_V0	  = 1/sigma*(diag(N)-1/sigma*H %*% R %*% solve(diag(q1)+1/sigma*t(H)%*%H%*%R)%*%t(H))
		inv_V0X 	= inv_V0%*%X                   #    To speed up the computation
		P 		  	= inv_V0-inv_V0X %*% solve(t(X)%*%inv_V0X)%*%t(inv_V0X)

    Tau_A_new = as.numeric(1/q1 * (Tau_A^2 * t(Y) %*% P %*% S1 %*% P %*% Y + Tau_A*q1 - sum(diag( Tau_A^2 * P %*% S1))))
    
    Y_star    = Y - Tau_A_new * S1 %*% P %*% Y
    
    IP_X      = diag(N)-P_X
    VA        =Tau_A*S1-(Tau_A*S1) %*% P %*% (Tau_A*S1)
    sigma_new = as.numeric(1/(N-p) * (t(Y_star) %*% IP_X %*% Y_star + sum(diag(IP_X%*%VA))))
    
    
		diff_A = abs((Tau_A_new-Tau_A)/Tau_A_new)		# use the relative difference rather than the absoult difference
		diff_s = abs((sigma_new-sigma)/sigma_new)
		if ((diff_A<0.0001)  & (diff_s<0.0001)) break
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
JointTest = function(Y,X,geno1,geno2){
	N 			= nrow(Y)															# the number of individual


	Simi1	 	= TypicalIBS(geno1)
# 	Simi1 = MaxIBS(geno1)
	S1			= Simi1$S

	Simi2	 	= TypicalIBS(geno2)
# 	Simi2 = MaxIBS(geno2)
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

#---------------------------- MAIN PART -------------------------------------  

setwd(path.expand("D:\\Projects\\Jung-Ying\\Real Data\\Alison Data"))

#---- Load data ----

Traits 		= as.matrix(read.table("Y.data",header=FALSE))
Y 				= array(Traits[,2], dim=c(301,1))		# Col 2 records BMI and Col 3 records log(BMI); it seems that using BMI would have a better result than log(BMI)
X					= as.matrix(read.table("X.data",header=FALSE))

geno1     = as.matrix(read.table("genoA.data",header=FALSE))
geno2     = as.matrix(read.table("genoB.data",header=FALSE))
# Jointp    = JointTest(Y,X,geno1,geno2)
# Epip      = EpiTest(Y,X,geno1,geno2)


  N 			= nrow(Y)															# the number of individual

	Simi1	= TypicalIBS(geno1)
	S1		= Simi1$S
	H1		= Simi1$H 
	R1		= Simi1$R
	q1		= Simi1$L

	Simi2	= TypicalIBS(geno2)
	S2		= Simi2$S
	H2		= Simi2$H 
	R2		= Simi2$R
	q2		= Simi2$L

	S12 	= S1 * S2														# The similarity matrix for the interaction

	p 		= ncol(X) 														# The number of fix effect

	P_X 	= X %*% solve(crossprod(X,X))%*% t(X) 				# The projection of X

	#----------------------- Initialize Tau and sigma ----------------------------

	Tau_A_new = 1                        
	Tau_B_new = 1
	sigma_new = 1
#   S2 = H2 %*% R2 %*% t(H2)
	#------------------------------- updata --------------------------------------

	repeat{		
		Tau_A		  = Tau_A_new
		Tau_B 	  = Tau_B_new
		sigma 	  = sigma_new
		
# 		R 			  = rbind(cbind(Tau_A*R1,array(0,dim = c(q1,q2))),cbind(array(0,dim = c(q2,q1)),Tau_B*R2))
# 		H 			  = cbind(H1,H2)
# 		V0			  = H %*% R %*% t(H) + sigma * diag(N)
    

# 		inv_V0	  = 1/sigma*(diag(N)-1/sigma*H %*% R %*% solve(diag(q1+q2)+1/sigma*t(H)%*%H%*%R)%*%t(H))
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
    
    cat("new ",Tau_A_new,Tau_B_new,sigma_new,"\n")
    
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






# Gene1p    = OneGeneTest(Y,X,geno2,geno1)
# Gene2p    = OneGeneTest(Y,X,geno1,geno2)

# cat("The Joint test for genoA and genoB is ",Jointp,"\n")
# cat("The Epi test for genoA and genoB is ",Epip,"\n")
# cat("The One gene test for gene1 is ", Gene1p, "\n")
# cat("The One gene test for gene2 is ", Gene2p, "\n")  