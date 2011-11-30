library(stats)
library(MASS)
library(CompQuadForm)
library(pls)
rm(list= ls())
set.seed(7)

#--------------------------- Load functions ----------------------------------

RealisticSimulation <- function(Sample_size,SNP_posi,Risk_Model, risk_effect){

  #------------------------------ 1.1. Load the data ------------------------------------

	Gene_RBJ 		    <- read.table("D:\\Projects\\Jung-Ying\\Codes\\Simulation\\type I error\\RBJ.txt", header=FALSE)      		#load the SNP information of the Gene RBJ
	Gene_RBJ			  <- as.matrix(Gene_RBJ)
	Gene_GPRC5B 	  <- read.table("D:\\Projects\\Jung-Ying\\Codes\\Simulation\\type I error\\GPRC5B.txt", header=FALSE)    	#load the SNP information of the Gene GPRC5B
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
		Main_effect <- 2*Causal_SNP[,1]
		Epi_effect  <- 0
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


LM_analysis <- function(Y,X,gene1,gene2){
  
    Inter = NULL
    
    for (i in 1:(ncol(gene1)-1)){
        for (j in i:(ncol(gene2))){
            Inter =  cbind(Inter,gene1[,i]*gene2[,j])
        }
    }   
    
    dat.all = data.frame(Y,X,gene1,gene2)
    fit0 = lm(Y~X,data = dat.all)     
    fit.all = update(fit0, . ~ . +gene2)
    p.all = anova(fit0,fit.all,test="F")
    return(p.all[2,6])
}
#---------------------------- Main part -------------------------------------

N     = 300
N.iter= 500


typeI = 0
#--- Generate all the possible SNP combinations for Model 1 --
SNP_posi = array(c(1,2,1,2),dim=c(1,4))


for (i in 1:N.iter){
  
  #- simulate the data ------
  SData = RealisticSimulation(Sample_size=N,SNP_posi,Risk_Model=0, risk_effect=0)
  Y     = SData$Y
  X     = SData$X
  gene1 = SData$geno1
  gene2 = SData$geno2
  
  
  #-- Using the LM--------
  pvalue = LM_analysis(Y,X,gene1,gene2)
  if (pvalue<0.05) typeI=typeI+1/N.iter
  cat(i,"\n")
}
    


