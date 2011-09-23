# Xin Wang, bioinformatics, NCSU, Aug, 2011

library(stats)
library(MASS)
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

weighted_chi2 <- function(eigenvalues,Ti){

	c1 	<- sum(eigenvalues)
	c2 	<- sum(eigenvalues^2)
	c3 	<- sum(eigenvalues^3)
	h 	<- c2^3/c3^2
	y 	<- (Ti-c1)*(h/c2)^(0.5)+h
	p 	<- pchisq(y,df=h)
	return(p_value=p)
	
	C 	<- rchisq(n=l*5000,df=1)
	dim(C) 	<- c(5000,l)	# so now the dimension of C is L*5000

	weightC <- C %*% Weights 
	N_w 		 <- nrow(weightC)
	p 			<- sum(as.numeric(T)<=weightC)/N_w
	
}

#-------------- Function : generate the similarity matrix from the genotype-----------
# This funcion is used to calculate the typical IBS from the given genotype

TypicalIBS <- function (geno){

	#- Sepetating the genotype into 2 haplotype matrix ----

	hA 	  		<- geno				# hA is used to how many A's in the allele, it equals to genotype
	ha 	  		<- 2-geno			# ha is used to how many a's in the allele, it equals to 2 - genotype
	
	N_gene  <- ncol(geno)    # the number of markers
	S_temp  <- hA %*% t(hA) + ha %*% t(ha)
	S         	<-  S_temp/(4*N_gene)
	
    eg        	<-  eigen(S, symmetric=T)
    evalue <-  eg$values
    le        	<-  evalue[1] # Biggest eigenvalue.
	if (le == 0){
		le 		<- 1
	}
	tmpkey <- (evalue / le) > 1e-7
	ll  	 		<- sum(tmpkey) # Rank of SSS
	RRRR	 	<- diag(evalue[tmpkey])
	HHHH   	<- eg$vectors[,tmpkey]
	return(list(S = S, R = RRRR, H = HHHH, L = ll))
	
}

ModifiedIBS <- function (geno){

	#- Sepetating the genotype into 2 haplotype matrix ----

	hA 	 		 <- geno				# hA is used to how many A's in the allele, it equals to genotype
	ha 	  		<- 2-geno			# ha is used to how many a's in the allele, it equals to 2 - genotype
	
	N 		  	<- nrow(geno)    # sample size
	N_gene <- ncol(geno)    # the number of markers
	
	S 		  	<- array(0,dim=c(N,N))
	
	#--- Generate the Sij -----
	
	for (i in 1:N){
		for(j in i:N){
			
			tempS 				<- (hA[i,]*hA[j,]+ha[i,]*ha[j,])/(4*N_gene)
			dim(tempS)  	<-  c(1,N_gene)
			S[i,j] 				<- sum(tempS)
			
			for (i1 in 1:(N_gene-1)){
				for (j1 in (i1+1):N_gene){
				
					S[i,j]		<- S[i,j]+tempS[1,i1]*tempS[1,j1]
					
				}
			}
			
			S[j,i] 				<- S[i,j]
		}
	}
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

EpiTest <- function(Y,genotype,N_gene1,N_gene2){

	N 			<- nrow(Y)															# the number of individual
	X 			<- array(1,dim=c(N,1))											# the population mean


	geno1 	<- genotype[,1:N_gene1]									# genotype for 1st gene
	geno2 	<- genotype[,(N_gene1+1):(N_gene1+N_gene2)]		# genotype for 2nd gene

	# Simi1	 	<- TypicalIBS(geno1)
	# Simi1 <- MaxIBS(geno1)
	Simi1    <-  ModifiedIBS(geno1)
	S1			<- Simi1$S
	H1			<- Simi1$H 
	R1			<- Simi1$R
	q1			<- Simi1$L

	# Simi2	 	<- TypicalIBS(geno2)
	# Simi2 <- MaxIBS(geno2)
	Simi2		<- ModifiedIBS(geno2) 
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
		inv_V0X 	<- inv_V0%*%X                   #    To speed up the computation
		P 		  	<- inv_V0-inv_V0X %*% solve(t(X)%*%inv_V0X)%*%t(inv_V0X)
		
		E_beta 		<- R %*% t(H) %*% P %*% Y
		E_betaA 	<- E_beta[1:q1]
		E_betaB 	<- E_beta[(q1+1):(q1+q2)]

		Var_beta 	<- R-R%*%t(H)%*%P%*%H%*%R
		Var_betaA <- Var_beta[1:q1,1:q1]
		Var_betaB <- Var_beta[(q1+1):(q1+q2),(q1+1):(q1+q2)] 

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
	V_eigen				<- e$vectors
	V_square_root	<- V_eigen %*% diag(sqrt(e$values)) %*% t(V_eigen)


	Weights_all		<- eigen(1/2*V_square_root %*% P0 %*% S12 %*% P0 %*% V_square_root, symmetric=TRUE, only.value=TRUE)
	temp					<- Weights_all$values
	temp2					<- sort(temp,decreasing=TRUE)
	dim(temp2)		<- c(N,1)
	big_enough 		<- sum(temp>10^-3)			# Get the number of big eigen values. here, the threshold for "big" is 10^-3 
	Weights				<- array(temp2[1:big_enough,1],dim=c(big_enough,1))

	C				     	<- rchisq(n=big_enough*5000,df=1)
	dim(C) 				<- c(5000,big_enough)	# so now the dimension of C is L*5000

	weightC				<- C %*% Weights 
	N_w				    <- nrow(weightC)
	p 						<- sum(as.numeric(T)<=weightC)/N_w
	return(P=p)
}

JointTest <- function(Y,genotype,N_gene1,N_gene2){
	N 			<- nrow(Y)															# the number of individual
	X 			<- array(1,dim=c(N,1))											# the population mean


	geno1 	<- genotype[,1:N_gene1]									# genotype for 1st gene
	geno2 	<- genotype[,(N_gene1+1):(N_gene1+N_gene2)]		# genotype for 2nd gene

	Simi1	 	<- TypicalIBS(geno1)
	#Simi1 <- MaxIBS(geno1)
	S1			<- Simi1$S

	Simi2	 	<- TypicalIBS(geno2)
	#Simi2 <- MaxIBS(geno2)
	S2			<- Simi2$S

	S12 		<- S1 * S2														# The similarity matrix for the interaction
	Q 		  <- diag(N) - X %*% solve(t(X) %*% X) %*% t(X)

	sigma 	<- as.numeric (t(Y) %*% Q %*% Y/(N-ncol(X)))

	P0 		  <- 1 / sigma * Q

	T 		  <- 1 / (2*sigma^2) * t(Y) %*% Q %*% (S1 + S2 + S12) %*% Q %*% Y
	Weights_all <- eigen(1 / (2 * sigma) * Q %*% (S1 + S2 + S12) %*% Q, symmetric=TRUE, only.value=TRUE)
	temp 	  <- Weights_all$values
	temp2 	<- sort(temp,decreasing=TRUE)
	dim(temp2) <- c(N,1)
	big_enough <- sum(temp>10^-3)			# Get the number of big eigen values. here, the threshold for "big" is 10^-3 

	Weights <- array(temp2[1:big_enough,1],dim=c(big_enough,1))


	C 		  <- rchisq(n=big_enough*5000,df=1)
	dim(C) 	<- c(5000,big_enough)	# so now the dimension of C is L*5000

	weightC <- C %*% Weights 
	N_w 	  <- nrow(weightC)
	p 		  <- sum(as.numeric(T)<=weightC)/N_w
	return(P=p)

}

#---------------------------- MAIN PART -------------------------------------  
	
#---- Load data ----

Traits 		<- as.matrix(read.table("good2_traits.txt",header=FALSE))
Y 					<- array(Traits[,4],dim=c(1068,1))
setwd(path.expand("D:/Written Prelim/Jung-Yng Zeng/data resort/Gene"))
AllGenes  	<-  list.files("D:/Written Prelim/Jung-Yng Zeng/data resort/Gene")
NoGene 	  	<-  length(AllGenes)

for (i in 1: (NoGene-1)){
	for (j in (i+1):NoGene){
	
		gene1 		<- as.matrix(read.table(AllGenes[i]))
		gene2 		<- as.matrix(read.table(AllGenes[j]))
		genotype<- cbind(gene1,gene2)
		N_gene1 	<- ncol(gene1)
		N_gene2 	<- ncol(gene2)
		Jointp 	<- JointTest(Y,genotype,N_gene1,N_gene2)
		Combi 		<- paste(AllGenes[i],AllGenes[j])
		Combi		<-  gsub(".txt", " ",Combi)
		cat(Combi,Jointp,"\n",file="result.txt",append=TRUE)
	}
}

# Epip		<- EpiTest(Y,genotype,N_gene1,N_gene2) 