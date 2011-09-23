# This Program is used to generate data to do joint test for our model
# We first try to simulate the data by haplotype bank, the simulation data would give u
# -s the phenotype, genotype, haplotype.

# Here, we try both the no weight version and weighted version for our method. The adva
# -ntage of weight is that it can help our method performs better under the situation t
# -hat the casual allele is a rara SNP than the no weighted version. However, its disad
# -vantage is that we lose some power when the causal one is a common SNP.

library(stats)
library(MASS)
library("pls")
library("CompQuadForm")
rm(list=ls())
set.seed(3)

#----------------------------------------------------------------------------------------------------------
#										FUNCTIONS														  |
#----------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------
#								1. SIMULATION										  |
#--------------------------------------------------------------------------------------


simulation 		<- function(Sample_size,SNP_posi,Risk_Model, risk_effect){

	#------------------------------ 1.1. Load the data ------------------------------------

	Gene_RBJ 		<- read.table("RBJ.txt", header=FALSE)      		#load the SNP information of the Gene RBJ
	Gene_RBJ			<- as.matrix(Gene_RBJ)
	Gene_GPRC5B 	<- read.table("GPRC5B.txt", header=FALSE)    	#load the SNP information of the Gene GPRC5B
	Gene_GPRC5B	<- as.matrix(Gene_GPRC5B)

	N_snp_RBJ  		<- ncol(Gene_RBJ)										#find how many SNPs in the RBJ gene
	N_snp_GPRC5B 	<- ncol(Gene_GPRC5B)									#find how many SNPs in the GPRC5B gene

	N_sample_RBJ	<- nrow(Gene_RBJ)										#find how many individuals in the RBJ gene
	N_sample_GPRC5B <- nrow(Gene_GPRC5B) 									#find how many individuals in the GPRC5B gene

	#------------------------------ 1.2. Set the parameters -------------------------------

	SNP11_posi	 	<- SNP_posi[1]                                    	    	#locate the 1st causal SNP in RBJ
	SNP12_posi		<- SNP_posi[2]												#locate the 2nd causal SNP in RBJ
	SNP21_posi	 	<- SNP_posi[3]                                    	    	#locate the 1st causal SNP in GPRC5B
	SNP22_posi		<- SNP_posi[4]												#locate the 2nd causal SNP in GPRC5B

	causal_posi 	<- c(SNP11_posi,SNP12_posi,SNP21_posi+N_snp_RBJ,SNP22_posi+N_snp_RBJ)		#when we consider the position as two genes, we need to add the number of first gene which is RBJ here

	#------------------------------ 1.3 Genearte the genotype -----------------------------

	Genotype 		<- array(0, dim=c(Sample_size, N_snp_RBJ+N_snp_GPRC5B))   	# the final genotype output
	Phenotype 		<- array(0, dim=c(Sample_size, 1))									# the final phenotype output
	
	tempA 			<- round(runif(Sample_size,1,N_sample_RBJ))					# randomly pick sample size individuals from the Gene bank of RBJ
	GeneA 			<- Gene_RBJ[tempA,]

	tempB 			<- round(runif(Sample_size,1,N_sample_GPRC5B))				# randomly pick sample size individuals from the Gene bank of GPRC5B
	GeneB 			<- Gene_GPRC5B[tempB,]

	genotype 		<- cbind(GeneA,GeneB)

	Causal_SNP		<- genotype[,causal_posi]
	
	if (Risk_Model==0){
		Main_effect <- 0
		Epi_effect 	 <-  0
	}
	
    if (Risk_Model==1){
		Main_effect <- risk_effect*(Causal_SNP[,1]+Causal_SNP[,3])
		Epi_effect  <- 0
	}
	
    if (Risk_Model==2){
		Main_effect <- 0
		Epi_effect  <- risk_effect*(Causal_SNP[,1]*Causal_SNP[,2]+Causal_SNP[,3]*Causal_SNP[,4])
	}
		
	if (Risk_Model==3){
		Main_effect <- risk_effect*Causal_SNP[,1]
		Epi_effect  <- risk_effect*(Causal_SNP[,2]*Causal_SNP[,4])
	}
	
	if (Risk_Model==4){
		Main_effect <- risk_effect*Causal_SNP[,1]
		Epi_effect  <- risk_effect*(Causal_SNP[,2]*Causal_SNP[,4]+Causal_SNP[,3]*Causal_SNP[,4])
	}
	
	if (Risk_Model==5){
		Main_effect <- 0
		Epi_effect  <- risk_effect*(Causal_SNP[,1]*Causal_SNP[,3]+Causal_SNP[,2]*Causal_SNP[,4])
	}

	#------------------------------ 1.4 Calculate the phenotype ---------------------------
	error_variance<- 1
	
	error		 		<- rnorm(Sample_size, 0,error_variance)

	Phenotype 		<- Main_effect + Epi_effect + error
	
	dim(Phenotype) <- c(Sample_size,1)

	Output 			 <- cbind(Phenotype, genotype)
	
	return(list(phenotype = Phenotype , genotype = genotype,N_gene1=N_snp_RBJ,N_gene2=N_snp_GPRC5B))
}

#--------------------------------------------------------------------------------------
#								2. Our methods										  |
#--------------------------------------------------------------------------------------
hsreg_VC <- function(Y,genotype,N_gene1,N_gene2){
	
		N 		<- nrow(Y)						# the number of individual
		X 		<- array(1,dim=c(N,1))			# the population mean
	
		#------------------ Using genotype to get the S -----------------------------------

		geno1 	<- genotype[,1:N_gene1]							# genotype for 1st gene
		geno2 	<- genotype[,(N_gene1+1):(N_gene1+N_gene2)]		# genotype for 2nd gene
		
		#- Sepetating the genotype into 2 haplotype matrix ----
		#-- for the geno1 --
		g1h1A_ <- array(0,dim=c(N,N_gene1))			# g1h1A_ is used to indicate whether the 1st allele is A for gene1
		g1h1_a <- array(0,dim=c(N,N_gene1))			# g1h1_a is used to indicate whether the 1st allele is a for gene1
		g1h2A_ <- array(0,dim=c(N,N_gene1))			# g1h2A_ is used to indicate whether the 2nd allele is A for gene1
		g1h2_a <- array(0,dim=c(N,N_gene1))			# g1h2_a is used to indicate whether the 2nd allele is a for gene1
		
		g1h1_a <- (geno1>0)*1						# This is a small algorithm to generate the g1h1A_~g1h2_a from the genotype
		g1h2_a <- geno1-g1h1_a

		g1h1A_ <- 1-g1h1_a
		g1h2A_ <- 1-g1h2_a
		
		#-- for the geno2 --
		g2h1A_ <- array(0,dim=c(N,N_gene2))			# This is similar to gene1
		g2h1_a <- array(0,dim=c(N,N_gene2))
		g2h2A_ <- array(0,dim=c(N,N_gene2))
		g2h2_a <- array(0,dim=c(N,N_gene2))
		
		g2h1_a <- (geno2>0)*1
		g2h2_a <- geno2-g2h1_a

		g2h1A_ <- 1-g2h1_a
		g2h2A_ <- 1-g2h2_a
		
		#----------------- Get the allele freq for each locus in each two gene-------------
		q1A 	<- array(1,dim=c(1,N_gene1))		# calculate the freq of main  allele for gene1
		q1a		<- array(1,dim=c(1,N_gene1))		# calculate the freq of minor allele for gene1
		q2A 	<- array(1,dim=c(1,N_gene2))		# calculate the freq of main  allele for gene2
		q2a		<- array(1,dim=c(1,N_gene2))		# calculate the freq of minor allele for gene2
		
		for (i in 1:N_gene1){
			q1A[i] <- (sum(geno1[,i]==0)*2+sum(geno1[,i]==1))/(2*N)
			q1a[i] <- 1-q1A[i]		
		}
		
		for (i in 1:N_gene2){
			q2A[i] <- (sum(geno2[,i]==0)*2+sum(geno2[,i]==1))/(2*N)
			q2a[i] <- 1-q2A[i]		
		}
		
		q1A[(q1A<=0.005)] <- 0.005
		q1A[(q1A>=0.995)] <- 0.995
		q1a[(q1a<=0.005)] <- 0.005
		q1a[(q1a>=0.995)] <- 0.995
		q2A[(q2A<=0.005)] <- 0.005
		q2A[(q2A>=0.995)] <- 0.995
		q2a[(q2a<=0.005)] <- 0.005
		q2a[(q2a>=0.995)] <- 0.995		
		
		w1A 	<-q1A^(-0.75)			# using the freq^(-3/4) as the weight
		w1a 	<-q1a^(-0.75)
		w2A 	<-q2A^(-0.75)
		w2a 	<-q2a^(-0.75)	
		
		S1_original <- array(0,dim=c(N,N))						# S by {0, 0.5, 1} value set
		S2_original <- array(0,dim=c(N,N))						# S by {0, 0.5, 1} value set
		
		# Following is small loop to generate the similarity matrix for gene1 and gene2	
		
		for (i in 1:N_gene1){
			temp1 <- w1A[i]*g1h1A_[,i] %*% t(g1h1A_[,i])+w1a[i]*g1h1_a[,i] %*% t(g1h1_a[,i]) + w1A[i]*g1h2A_[,i] %*% t(g1h2A_[,i])+w1a[i]*g1h2_a[,i] %*% t(g1h2_a[,i])
			temp2 <- w1A[i]*g1h1A_[,i] %*% t(g1h2A_[,i])+w1a[i]*g1h1_a[,i] %*% t(g1h2_a[,i]) + w1A[i]*g1h2A_[,i] %*% t(g1h1A_[,i])+w1a[i]*g1h2_a[,i] %*% t(g1h1_a[,i])
			S1_original <- S1_original+ temp1*(temp1>=temp2)+temp2*(temp1<temp2)
		}
		
		for (i in 1:N_gene2){
			temp1 <- w2A[i]*g2h1A_[,i] %*% t(g2h1A_[,i])+w2a[i]*g2h1_a[,i] %*% t(g2h1_a[,i]) + w2A[i]*g2h2A_[,i] %*% t(g2h2A_[,i])+w2a[i]*g2h2_a[,i] %*% t(g2h2_a[,i])
			temp2 <- w2A[i]*g2h1A_[,i] %*% t(g2h2A_[,i])+w2a[i]*g2h1_a[,i] %*% t(g2h2_a[,i]) + w2A[i]*g2h2A_[,i] %*% t(g2h1A_[,i])+w2a[i]*g2h2_a[,i] %*% t(g2h1_a[,i])
			S2_original <- S2_original+ temp1*(temp1>=temp2)+temp2*(temp1<temp2)
		}
		
		S1 <- S1_original/(2*N_gene1)			
		S2 <- S2_original/(2*N_gene2)
		S12 <- S1* S2
		
		Q 		<- diag(N) - X %*% solve(t(X) %*% X) %*% t(X)
		
		sigma 	<- as.numeric (t(Y) %*% Q %*% Y/(N-ncol(X)))

		P0 		<- 1 / sigma * Q
		
		T 		<- 1 / (2*sigma^2) * t(Y) %*% Q %*% (S1 + S2 + S12) %*% Q %*% Y
		Weights_all <- eigen(1 / (2 * sigma) * Q %*% (S1 + S2 + S12) %*% Q, symmetric=TRUE, only.value=TRUE)
		temp 	<- Weights_all$values
		temp2 	<- sort(temp,decreasing=TRUE)
		dim(temp2) <- c(N,1)
		big_enough <- sum(temp>10^-3)			# Get the number of big eigen values. here, the threshold for "big" is 10^-3 
		
		Weights <- array(temp2[1:big_enough,1],dim=c(big_enough,1))
		
		
		C 		<- rchisq(n=big_enough*5000,df=1)
		dim(C) 	<- c(5000,big_enough)	# so now the dimension of C is L*5000

		weightC <- C %*% Weights 
		N_w 	<- nrow(weightC)
		p 		<- sum(as.numeric(T)<=weightC)/N_w
	
	# p <- farebrother(as.numeric(T),Weights)$res

    return(p_value=p)
}

#--------------------------------------------------------------------------------------
#                           3. Use PCA to fit the model                               |
#--------------------------------------------------------------------------------------

 PCA_analysis <- function(phenotype,genotype,n_gene1,n_gene2){

	N 		<- nrow(phenotype)								# the nmuber of individuals in the sample
	
	geno1 	<- genotype[,1:n_gene1]							# the genotype for the 1st gene
	geno2 	<- genotype[,(n_gene1+1):(n_gene1+n_gene2)]		# the genotype for the 2nd gene
	
	#---------- 2.1 get the 1st PCA component of geno1 and geno2)----------------------

	gene1_PCA		 <- princomp(geno1,cor=FALSE, scores=TRUE)
	gene2_PCA		 <- princomp(geno2,cor=FALSE, scores=TRUE)
	
	Z1				 <- gene1_PCA$score[,1]
	Z2				 <- gene2_PCA$score[,1]
	
	model_analysis 	 <- glm(phenotype~genotype+Z1*Z2)		# the alternative model
	
	ptest1		 	<- anova(model_analysis,test="Chisq")
	po1 		 	<- 1-pchisq((max(ptest1[[4]])-min(ptest1[[4]])),df=max(ptest1[[3]])-min(ptest1[[3]]))
	return(po1)
}

#--------------------------------------------------------------------------------------
#                           4. Use PLS to fit the model                               |
#--------------------------------------------------------------------------------------

 PLS_analysis <- function(phenotype,genotype,n_gene1,n_gene2){

	N 		<- nrow(phenotype)								# the nmuber of individuals in the sample
	
	geno1 	<- genotype[,1:n_gene1]							# the genotype for the 1st gene
	geno2 	<- genotype[,(n_gene1+1):(n_gene1+n_gene2)]		# the genotype for the 2nd gene
	
	#---------- 2.1 get the 1st PLS component regression geno1 on (phenotype, geno2)---
	
	center_geno1 <- apply(geno1, 2, function(x)(x-mean(x))/sd(x))
	center_geno2 <- apply(geno2, 2, function(x)(x-mean(x))/sd(x))
	
	#------------------------ The following 4 lines are used to remove NaN ------------
	bad			 <- sapply(center_geno1[1,], function(x) all(is.nan(x)))
	center_geno1 <- center_geno1[,!bad]
	bad			 <- sapply(center_geno2[1,], function(x) all(is.nan(x)))
	center_geno2 <- center_geno2[,!bad]	

	Y 			 <- (phenotype-mean(phenotype))/sd(phenotype)	
	
	dat1		 <- data.frame(phenotype=phenotype,center_geno2)
	fit1		 <- glm(phenotype~.,data=dat1)
	mu			 <- fitted.values(fit1)
	
	dat2 		 <- data.frame(Gene1=I(center_geno1),Y_Gene2=I(cbind(center_geno2,Y)))
	pls1		 <- plsr(Y_Gene2 ~ Gene1, data = dat2)
	In			 <- scores(pls1)[,1]*mu
	dat3		 <- data.frame(phenotype=phenotype,geno1,geno2,In)
	fit3		 <- glm(phenotype~.,data=dat3)
	ptest1		 <- anova(fit3,test="Chisq")
	po1 		 <- 1-pchisq((max(ptest1[[4]])-min(ptest1[[4]])),df=max(ptest1[[3]])-min(ptest1[[3]]))
	
	return(po1)
}

#----------------------------------------------------------------------------------------------------------
#										Main part														  |
#----------------------------------------------------------------------------------------------------------

N 			<- 800

#--------------------------------- Load data by different effect ------------------------------------------
SNP1 <- array(c(1,2,3),dim=c(1,3))
SNP2 <- array(c(5,6,8),dim=c(1,3))
SNP_posi=array(0,dim=c(36,4))

p <- 1
for (i1 in 1:3 ){
	
	for (i2 in SNP1[-i1]){
		for (j1 in 1:3){
			temp2 <- 1:3
			for (j2 in temp2[-j1]){
				SNP_posi[p,] <- array(c(SNP1[i1],SNP1[i2],SNP2[j1],SNP2[j2]),dim=c(1,4))
				p <- p+1
			}
		}
	}
}

Risk_Model <- 2
risk_effect <- 0.15
iteration 	<- 200

for (i in 1:36){
	P_value		<- array(0,dim=c(3,iteration))

	for (iter in 1:iteration){

		# ------------------------------- Generate the data --------------------------------

		all_data <- simulation(N,SNP_posi[i,],Risk_Model, risk_effect)
		phenotype<- all_data$phenotype
		genotype <- all_data$genotype
		N_gene1  <- all_data$N_gene1			# number of SNPs in the 1st gene
		N_gene2  <- all_data$N_gene2			# number of SNPs in the 2nd gene

		# -------------------------------- Analysis the data -------------------------------

		#P_value[1, iter] <- hsreg_VC(phenotype,genotype,N_gene1,N_gene2)
		P_value[2, iter] <- PCA_analysis(phenotype,genotype,N_gene1,N_gene2)
		P_value[3, iter] <- PLS_analysis(phenotype,genotype,N_gene1,N_gene2)
	}
	p_val <- rowSums(P_value<=0.05)/iteration
	cat("SNP_posi",SNP_posi[i,]," Risk Model ",Risk_Model, " p_value ", p_val[1],p_val[2],p_val[3],"\n")
}