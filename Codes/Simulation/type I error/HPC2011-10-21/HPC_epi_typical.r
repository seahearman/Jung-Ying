library(stats)
library(MASS)
library(CompQuadForm)
rm(list= ls())
set.seed(4)

#--------------------------- Load functions ----------------------------------

RealisticSimulation <- function(Sample_size,SNP_posi,Risk_Model, risk_effect){

  #------------------------------ 1.1. Load the data ------------------------------------

	Gene_RBJ 		    <- read.table("RBJ.txt", header=FALSE)      		#load the SNP information of the Gene RBJ
	Gene_RBJ			  <- as.matrix(Gene_RBJ)
	Gene_GPRC5B 	  <- read.table("GPRC5B.txt", header=FALSE)    	#load the SNP information of the Gene GPRC5B
	Gene_GPRC5B	    <- as.matrix(Gene_GPRC5B)

	N_snp_RBJ  		  <- ncol(Gene_RBJ)										#find how many SNPs in the RBJ gene
	N_snp_GPRC5B 	  <- ncol(Gene_GPRC5B)									#find how many SNPs in the GPRC5B gene

	N_sample_RBJ	  <- nrow(Gene_RBJ)										#find how many individuals in the RBJ gene
	N_sample_GPRC5B <- nrow(Gene_GPRC5B) 									#find how many individuals in the GPRC5B gene

	#------------------------------ 1.2. Set the parameters -------------------------------

	SNP11_posi	 	<- SNP_posi[1]                                    	    	#locate the 1st causal SNP in RBJ
	SNP12_posi		<- SNP_posi[2]												#locate the 2nd causal SNP in RBJ
	SNP21_posi	 	<- SNP_posi[3]                                    	    	#locate the 1st causal SNP in GPRC5B
	SNP22_posi		<- SNP_posi[4]												#locate the 2nd causal SNP in GPRC5B

	causal_posi 	<- c(SNP11_posi,SNP12_posi,SNP21_posi+N_snp_RBJ,SNP22_posi+N_snp_RBJ)		#when we consider the position as two genes, we need to add the number of first gene which is RBJ here

	#------------------------------ 1.3 Genearte the genotype -----------------------------

	Genotype 		<- array(0, dim=c(Sample_size, N_snp_RBJ+N_snp_GPRC5B))   	# the final genotype output
	Phenotype 	<- array(0, dim=c(Sample_size, 1))									# the final phenotype output
	
	tempA 			<- round(runif(Sample_size,1,N_sample_RBJ))					# randomly pick sample size individuals from the Gene bank of RBJ
	GeneA 			<- Gene_RBJ[tempA,]

	tempB 			<- round(runif(Sample_size,1,N_sample_GPRC5B))				# randomly pick sample size individuals from the Gene bank of GPRC5B
	GeneB 			<- Gene_GPRC5B[tempB,]

	genotype 		<- cbind(GeneA,GeneB)

	Causal_SNP	<- genotype[,causal_posi]
	
	if (Risk_Model==0){
		Main_effect <- 0
		Epi_effect 	<- 0
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
  
  X = array(1,dim=c(Sample_size,1))
  
  geno1        <- genotype[,1:N_snp_RBJ]
  geno2        <- genotype[,-(1:N_snp_RBJ)]
	
	return(list(Y = Phenotype , X=X,geno1 = geno1,geno2 = geno2))
}
TypicalIBS_genotype     = function (geno){
  
  #- If the data has a 2-allelic genotype format, i.e. there are only 2 allele for each
  #- marker, say A, a, and the genotype is recorded by a scale, that is, AA:0, Aa:1 and aa:2

  #- if the genotype has more than 2 allelic data, we may want to calculated in other way
  #- Sepetating the genotype into 2 haplotype matrix ----
  #---------------------------------------------------------
  # this is the version can be applied to the 2-allelic data
  #---------------------------------------------------------
  N_gene = ncol(geno)
	N 		 = nrow(geno)

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

EpiTest = function(Y,X,Simi1,Simi2){

  N       = nrow(Y)								# the number of individual

	S1		= Simi1$S
	q1		= Simi1$L

	S2		= Simi2$S
	q2		= Simi2$L

	S12 	= S1 * S2									# The similarity matrix for the interaction

	p 		= ncol(X) 								# The number of fix effect

	P_X 	= X %*% solve(crossprod(X,X))%*% t(X) 				# The projection of X

	#----------------------- Initialize Tau and sigma ----------------------------

	Tau_A_new = 0.5                        
	Tau_B_new = 0.5
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
    
  	diff_A = abs((Tau_A_new-Tau_A)/Tau_A_new)		# use the relative difference rather than the absoult difference
		diff_B = abs((Tau_B_new-Tau_B)/Tau_B_new)
		diff_s = abs((sigma_new-sigma)/sigma_new)

		if ((diff_A<0.001) & (diff_B<0.001) & (diff_s<0.001)) break
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
#---------------------------- Main part -------------------------------------

N     = 300
N.iter= 500
SNP_posi=array(c(1,2,1,2),dim=c(1,4))
p     = 0
for (i in 1:N.iter){
  SData = RealisticSimulation(Sample_size=N,SNP_posi,Risk_Model=0, risk_effect=0)
  Y     = SData$Y
  X     = SData$X
  gene1 = SData$geno1
  gene2 = SData$geno2
  
  Simi1     = TypicalIBS_genotype(gene1)
  Simi2     = TypicalIBS_genotype(gene2)
  
  pvalue    = EpiTest(Y,X,Simi2,Simi1)
  cat("round=",i,"\n")
  if (pvalue<0.05) p=p+1
}
type.I.error = p/N.iter
cat("The type I error for Epi test is ",type.I.error,"\n")