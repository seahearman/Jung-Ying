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
