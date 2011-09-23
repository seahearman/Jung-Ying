# Xin Wang, bioinformatics, NCSU, Aug, 2011

library(stats)
library(MASS)
rm(list= ls())
set.seed(4)

#--------------------------------------------------------------------------------------------------
#														FUNCTIONS																		|
#--------------------------------------------------------------------------------------------------

#-------------------------FUNCTION 1: SIMULATION -------------------------------
 
# The simulation function is to genearate data accroding to the model:
# Y=X*gammar+H^A * beta^A + H^B * beta^B + e

simulation <- function(N, Tau_A,Tau_B,phi, sigma){
	
	#1. Generate the fix effect
	
	X 		<- array(1, dim <- c(N,2))
	X[,2] <- rnorm(N,mean=0, sd=1)
	Gamma <- c(1,1)

	#2 Generate the random effect
	#2.1. Generate the H1 and H2 : haplotype matrix
	
	p 		<- c(0.397,0.208,0.02,0.013,0.018,0.138,0.129,0.054)  # the propability of the 8 haplotypes in population
	
	#---- Generate the haplotype for Gene A
	
	H1 	<- t(rmultinom(N,2,p))
	
	#---- Generate the haplotype for Gene B

	H2 	<- t(rmultinom(N,2,p))
	
	#---- Generate the genotype from the haplotype
	
	H 		<- array(0,dim <- c(8,5))
	H[1,] <- c(0, 0, 0, 0, 0)
	H[2,] <- c(0, 0, 0, 0, 1)
	H[3,] <- c(0, 1, 0, 0, 0)
	H[4,] <- c(0, 1, 0, 0, 1)
	H[5,] <- c(0, 1, 1, 0, 0)
	H[6,] <- c(0, 1, 1, 1, 1)
	H[7,] <- c(1, 0, 0, 0, 0)
	H[8,] <- c(1, 0, 0, 0, 1)
	
	SA	    <- H1 %*% H
	SB 	<- H2 %*% H	
	
	qA 	<- length(p)
	qB 	<- length(p)
	
	H_AB 	<- array(0, dim <- c(N,qA*qB))
	# 2.2 Genearte the H_AB from H1 and H2
	pcol 	<- 0
	for (A in 1:qA){
	
		for (B in 1:qB){
		
			pcol 			 <- pcol+1
			H_AB[, pcol] <- H1[,A]*H2[,B]
		
		}
	
	}

	#2.3 Generate the variance of beta: R1 and R2    

	R_all <- c( 1,0.8,0.8,0.6,0.6,0.2,0.8,0.6,
					0.8,1,0.6,0.8,0.4,0.4,0.6,0.8,
					0.8,0.6,1,0.8,0.8,0.4,0.6,0.4,
					0.6,0.8,0.8,1,0.6,0.6,0.4,0.6,
					0.6,0.4,0.8,0.6,1,0.6,0.4,0.2,
					0.2,0.4,0.4,0.6,0.6,1.0,0,0.2,
					0.8,0.6,0.6,0.4,0.4,0,1.0,0.8,
					0.6,0.8,0.4,0.6,0.2,0.2,0.8,1)   # the similarity matrix of the 8 haplotypes

	dim(R_all) <- c(8,8)

	R1 <- R_all
	R2 <- R_all
	
	#2.4 generate R_AB according to R1 and R2 (there are 4 loops, really ineffiecent)
	
	R_AB <- array(0, dim <- c(qA*qB,qA*qB))
	prow <- 0
	for (A1 in 1:qA){
	
		for (B1 in 1:qB){
		
			prow <- prow+1
			pcol <- 0
			for (A2 in 1:qA){
			
				for (B2 in 1:qB){
				
					pcol 				   <- pcol+1
					R_AB[prow, pcol] <- R1[A1,A2]*R2[B1,B2]
					
				}
				
			}
			
		}
		
	}
	Sigma_beta_A 	<- Tau_A*R1
	Sigma_beta_B 	<- Tau_B*R2
	Sigma_beta_AB 	<- phi*R_AB

	beta_A 			<- mvrnorm(n = 1, mu=rep(0,qA), Sigma_beta_A, tol = 1e-6, empirical = FALSE)
	beta_B 			<- mvrnorm(n = 1, mu=rep(0,qB), Sigma_beta_B, tol = 1e-6, empirical = FALSE)
	beta_AB 			<- mvrnorm(n = 1, mu=rep(0,qA*qB), Sigma_beta_AB, tol = 1e-6, empirical = FALSE)
	
	#3 Generate the error term
	e 					<- rnorm(N,mean=0,sd=sqrt(sigma))

	#4. sum everything up!	
	Y 					<- X %*% Gamma+ H1 %*% beta_A + H2%*%beta_B+H_AB %*% beta_AB+e
	
	return(list(Y = Y, X = X, genotype = cbind(SA,SB), N_gene1 = 5, N_gene2 = 5))
}

#-------------------- Function 2: weighted chi-square distribution---------------------
# The main idea for for this function is:
# 1. generate 8*10,000 chi square distributed random variables: Matrix C(this should be
#    outside of the function to save time, you don't need to generate lots of random va
#    -irables each time)
# 2. Acoording to the eigenvalues, calculated 10,000 weighted chi-square distributed ra
#    -ndom variables
# 3. given the Ti statistics, get the p-value.

weighted_chi2 <- function(eigenvalues,Ti){

	c1 <- sum(eigenvalues)
	c2 <- sum(eigenvalues^2)
	c3 <- sum(eigenvalues^3)
	h 	<- c2^3/c3^2
	y 	<- (Ti-c1)*(h/c2)^(0.5)+h
	p 	<- pchisq(y,df=h)
	return(p_value=p)
	
	l <- 
	C 		<- rchisq(n=l*5000,df=1)
	dim(C) 	<- c(5000,l)	# so now the dimension of C is L*5000

	weightC <- C %*% Weights 
	N_w 	<- nrow(weightC)
	p 		<- sum(as.numeric(T)<=weightC)/N_w
	
}

#-------------- Function 3: generate the similarity matrix from the genotype-----------
# This funcion is used to calculate the typical IBS from the given genotype

TypicalIBS <- function (geno){

	#- Sepetating the genotype into 2 haplotype matrix ----

	hA 	  <- geno				# hA is used to how many A's in the allele, it equals to genotype
	ha 	  <- 2-geno			# ha is used to how many a's in the allele, it equals to 2 - genotype
	
	S_temp  <- hA %*% t(hA) + ha %*% t(ha)
	S         <-  S_temp/4
	
    eg        <-  eigen(S, symmetric=T)
    evalue  <-  eg$values
    le        <-  evalue[1] # Biggest eigenvalue.
	if (le == 0){
		le <- 1
	}
	tmpkey <- (evalue / le) > 1e-7
	ll  	 <- sum(tmpkey) # Rank of SSS
	RRRR	 <- diag(evalue[tmpkey])
	HHHH   <- eg$vectors[,tmpkey]
	return(list(S = S, R = RRRR, H = HHHH, L = ll))
	
}

MaxIBS		 <- function (geno){

	N_gene <- ncol(geno)
	
	h1A_ <- array(0,dim=c(N,N_gene))			# h1A_ is used to indicate whether the 1st allele is A for gene
	h1_a <- array(0,dim=c(N,N_gene))			# h1_a is used to indicate whether the 1st allele is a for gene
	h2A_ <- array(0,dim=c(N,N_gene))			# h2A_ is used to indicate whether the 2nd allele is A for gene
	h2_a <- array(0,dim=c(N,N_gene))			# h2_a is used to indicate whether the 2nd allele is a for gene

	h1_a <- (geno>0)*1						# This is a small algorithm to generate the h1A_~h2_a from the genotype
	h2_a <- geno-h1_a

	h1A_ <- 1-h1_a
	h2A_ <- 1-h2_a


	#----------------- Get the allele freq for each locus in each two gene-------------
	qA 	<- array(1,dim=c(1,N_gene))		# calculate the freq of main  allele for gene
	qa		<- array(1,dim=c(1,N_gene))		# calculate the freq of minor allele for gene



	wA 	<-qA^(-0.75)			# using the freq^(-3/4) as the weight
	wa 	<-qa^(-0.75)

	S_original <- array(0,dim=c(N,N))						# S by {0, 0.5, 1} value set


	# Following is small loop to generate the similarity matrix for gene1 and gene2	

	for (i in 1:N_gene){
		temp1 <- wA[i]*h1A_[,i] %*% t(h1A_[,i])+wa[i]*h1_a[,i] %*% t(h1_a[,i]) + wA[i]*h2A_[,i] %*% t(h2A_[,i])+wa[i]*h2_a[,i] %*% t(h2_a[,i])
		temp2 <- wA[i]*h1A_[,i] %*% t(h2A_[,i])+wa[i]*h1_a[,i] %*% t(h2_a[,i]) + wA[i]*h2A_[,i] %*% t(h1A_[,i])+wa[i]*h2_a[,i] %*% t(h1_a[,i])
		S_original <- S_original+ temp1*(temp1>=temp2)+temp2*(temp1<temp2)
	}

	S <- S_original/(2*N_gene)
	
    eg        <-  eigen(S, symmetric=T)
    evalue  <-  eg$values
    le        <-  evalue[1] # Biggest eigenvalue.
	if (le == 0){
		le <- 1
	}
	tmpkey <- (evalue / le) > 1e-7
	ll  	 <- sum(tmpkey) # Rank of SSS
	RRRR	 <- diag(evalue[tmpkey])
	HHHH   <- eg$vectors[,tmpkey]
	return(list(S = S, R = RRRR, H = HHHH, L = ll))	

}

#---------------------------- MAIN PART -------------------------------------  
	
#---- Load data ----
Traits 	<- as.matrix(read.table("good2_traits.txt",header=FALSE))
Y 			<- array(Traits[,4],dim=c(1068,1))
genotype <- as.matrix(read.table("gain2.txt",header=FALSE))
N_gene1 	<- 2
N_gene2 	<- 2

N 			<- nrow(Y)															# the number of individual
X 			<- array(1,dim=c(N,1))											# the population mean


geno1 	<- genotype[,1:N_gene1]									# genotype for 1st gene
geno2 	<- genotype[,(N_gene1+1):(N_gene1+N_gene2)]		# genotype for 2nd gene

Simi1	 	<- TypicalIBS(geno1)
#Simi1 <- MaxIBS(geno1)
S1			<- Simi1$S
H1			<- Simi1$H 
R1			<- Simi1$R
q1			<- Simi1$L

Simi2	 	<- TypicalIBS(geno2)
#Simi2 <- MaxIBS(geno2)
S2			<- Simi2$S
H2			<- Simi2$H 
R2			<- Simi2$R
q2			<- Simi2$L

S12 		<- S1 * S2														# The similarity matrix for the interaction

p 			<- ncol(X) 														# The number of fix effect
N			<- nrow(X) 														# Sample size

P_X 		<- X %*% solve(crossprod(X,X))%*% t(X) 				# The projection of X

#----------------------- Initialize Tau and sigma ----------------------------

Tau_A_new <- 1                        
Tau_B_new <- 1
sigma_new <- 1

#------------------------------- updata --------------------------------------

repeat{		
	Tau_A		  <- as.numeric(Tau_A_new)
	Tau_B 	  <- as.numeric(Tau_B_new)
	sigma 	  <- as.numeric(sigma_new)
	
	R 			  <- rbind(cbind(Tau_A*R1,array(0,dim <- c(q1,q2))),cbind(array(0,dim <- c(q2,q1)),Tau_B*R2))
	H 			  <- cbind(H1,H2)
	V0			  <- H %*% R %*% t(H) + sigma * diag(N)

	inv_V0	  <- 1/sigma*(diag(N)-1/sigma*H %*% R %*% solve(diag(q1+q2)+1/sigma*t(H)%*%H%*%R)%*%t(H))
	inv_V0X 	  <- inv_V0%*%X                   #    To speed up the computation
	P 		  	  <- inv_V0-inv_V0X %*% solve(t(X)%*%inv_V0X)%*%t(inv_V0X)
	
	E_beta 			<- R %*% t(H) %*% P %*% Y
	E_betaA 			<- E_beta[1:q1]
	E_betaB 			<- E_beta[(q1+1):(q1+q2)]

	Var_beta 		<- R-R%*%t(H)%*%P%*%H%*%R
	Var_betaA 		<- Var_beta[1:q1,1:q1]
	Var_betaB 		<- Var_beta[(q1+1):(q1+q2),(q1+1):(q1+q2)] 

	inv_RE_betaA 	<- Tau_A*t(H1)%*%P%*%Y
	inv_RVar_betaA<- Tau_A*diag(q1)-Tau_A^2*t(H1)%*%P%*%H1%*%R1

	inv_RE_betaB 	<- Tau_B*t(H2)%*%P%*%Y
	inv_RVar_betaB<- Tau_B*diag(q2)-Tau_B^2*t(H2)%*%P%*%H2%*%R2

	Tau_A_new 		<- 1/q1*(t(E_betaA)%*%inv_RE_betaA+sum(diag(inv_RVar_betaA)))
	Tau_B_new 		<- 1/q2*(t(E_betaB)%*%inv_RE_betaB+sum(diag(inv_RVar_betaB)))
	sigma_new 		<- 1/(N-p)*(crossprod((Y-H%*%E_beta),(diag(N)-P_X))%*% (Y-H%*%E_beta)+sum(diag(crossprod(H, (diag(N)-P_X))%*%H%*%Var_beta)))

	diff_A <- abs((as.numeric(Tau_A_new)-Tau_A)/as.numeric(Tau_A_new))		# use the relative difference rather than the absoult difference
	diff_B <- abs((as.numeric(Tau_B_new)-Tau_B)/as.numeric(Tau_B_new))
	diff_s <- abs((as.numeric(sigma_new)-sigma)/as.numeric(sigma_new))
	if ((diff_A<0.001) & (diff_B<0.001) & (diff_s<0.001)) break
}

P0 					  <- P 

T0 					  <- 1/2*t(Y) %*% P0 %*% S12 %*% P0 %*% Y	#Get T0 under null hypothesis
e						  <- eigen(V0, symmetric=TRUE)
V_eigen				  <- e$vectors
V_square_root		  <- V_eigen %*% diag(sqrt(e$values)) %*% t(V_eigen)


Weights_all		  <- eigen(1/2*V_square_root %*% P0 %*% S12 %*% P0 %*% V_square_root, symmetric=TRUE, only.value=TRUE)
temp					  <- Weights_all$values
temp2					  <- sort(temp,decreasing=TRUE)
dim(temp2)			  <- c(N,1)
big_enough 			  <- sum(temp>10^-3)			# Get the number of big eigen values. here, the threshold for "big" is 10^-3 
Weights				  <- array(temp2[1:big_enough,1],dim=c(big_enough,1))

C				     	  <- rchisq(n=big_enough*5000,df=1)
dim(C) 				  <- c(5000,big_enough)	# so now the dimension of C is L*5000

weightC				  <- C %*% Weights 
N_w				      <- nrow(weightC)
p 						  <- sum(as.numeric(T)<=weightC)/N_w
