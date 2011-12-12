library(stats)
library(MASS)
library(CompQuadForm)
library(pls)
rm(list= ls())
set.seed(7)

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
  	
	
	if (Risk_Model==1){
		Main_effect <- 0
		Epi_effect  <- risk_effect*(Causal_SNP[,1]*Causal_SNP[,3]+Causal_SNP[,2]*Causal_SNP[,4])
	}	

  if (Risk_Model==2){
		Main_effect <- risk_effect*(Causal_SNP[,1]+Causal_SNP[,3])
		Epi_effect  <- 0
	}
	
    if (Risk_Model==3){
		Main_effect <- 0
		Epi_effect  <- risk_effect*(Causal_SNP[,1]*Causal_SNP[,2]+Causal_SNP[,3]*Causal_SNP[,4])
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

JY.Typical.IBS = function(SNP=addsnp, weight.power = (-1), Gmode="add", MAF=NULL)

{

  ## tSNP is SNPs-rows, individuals=columns format.
  
  tSNP = t(SNP)
  
  nsubj = ncol(tSNP)
  
  nloci = nrow(tSNP)
  
  if(Gmode=="add")
  
  {
  
    ## Weights for the S matrix calculated from the allele frequencies.
    
    tSNP.A = tSNP
    
    tSNP.a = 2 - tSNP.A
    
    AFA = rowSums(tSNP.A) / (2 * nsubj)
    
    AFa = 1 - AFA # Allele Freq of "a"
    
    wt.A = AFA ^ weight.power;
    
    wt.a = AFa ^ weight.power;
    
    ## use the following if use Wu weights
    
    ## wt.a = wt.A
    
    wttSNP.A = sqrt(wt.A) * tSNP.A
    
    wttSNP.a = sqrt(wt.a) * tSNP.a
    
    ## This is the avgIBS weighted S matrix.
    
    SS.wt = ( crossprod(wttSNP.A) + crossprod(wttSNP.a) ) / nloci
    
    ##-------------------
    
    ## SSmat by maxIBS kernel
    
    ##-------------------
    
    ## The maxIBS values can be calculated from the avgIBS S matrix
    
    ## plus an extra term for the het-het cases.
    
    mg2 = 1 - (1-tSNP)^2
    
    wtmg2.A = sqrt(wt.A) * mg2
    
    wtmg2.a = sqrt(wt.a) * mg2
    
    SS.hh = (crossprod(wtmg2.A) + crossprod(wtmg2.a)) / nloci
    
    SS.maxibs = SS.wt + SS.hh
    
    S      = SS.maxibs
    eg     =  eigen(S, symmetric=T)
    evalue =  eg$values
    le     =  evalue[1] # Biggest eigenvalue.
    if (le == 0){
      le     = 1
    }
    tmpkey = (evalue / le) > 1e-7
  	l.maxibs = sum(tmpkey) # Rank of SSS
  
  }
  
  return(list("S"=SS.maxibs,"L"=l.maxibs))

}
JointTest = function(Y,X,Simi1,Simi2){
  N   		= nrow(Y)															# the number of individual

	S1			= Simi1$S

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

#---------------------------- Main part -------------------------------------

N     = 300
N.iter= 100

GeneA  = array(c(1,3,7),dim=c(1,3))
GeneA  = array(c(1,3,7),dim=c(1,3))
ALD    = array(c(0.588,0.126,0.159),dim=c(1,3))
AMAF   = array(c(0.115,0.482,0.062),dim=c(1,3))
LableA = array(c("HM","LC","LR"),dim=c(1,3))
GeneB  = array(c(7,9,1,6,3),dim=c(1,5))
BLD    = array(c(0.259,0.23,0.186,0.136,0.065),dim=c(1,5))
BMAF   = array(c(0.46,0.482,0.146,0.159,0.044),dim=c(1,5))
LableB = array(c("HC","NC","NM","LM","LR"),dim=c(1,5))

#--- Generate all the possible SNP combinations for Model 1 --
SNP_posi  = NULL
Lable     = NULL
LDPatten  = NULL
MAFPatten = NULL

for (gA1 in 1:2){
  
  SA1   = GeneA[1,gA1]
  LDA1  = ALD[1,gA1]
  MAFA1 = AMAF[1,gA1]
  LA1   = LableA[1,gA1]
  
  for (gA2 in (gA1+1):3){
    
    SA2   = GeneA[1,gA2]
    LA2   = LableA[1,gA2]
    LDA2  = ALD[1,gA2]
    MAFA2 = AMAF[1,gA2]
    
    for (gB1 in 1:4){
      
      SB1   = GeneB[1,gB1]
      LB1   = LableB[1,gB1]
      LDB1  = BLD[1,gB1]
      MAFB1 = BMAF[1,gB1]
      
      for (gB2 in (gB1+1):5){
        
        SB2   = GeneB[1,gB2]
        LB2   = LableB[1,gB2]
        LDB2  = BLD[1,gB2]
        MAFB2 = BMAF[1,gB2]
        
        SNP_posi  = rbind(SNP_posi,array(c(SA1,SA2,SB1,SB2),dim=c(1,4)))
        LDPatten  = rbind(LDPatten,array(c(LDA1,LDA2,LDB1,LDB2),dim=c(1,4)))
        MAFPatten = rbind(MAFPatten,array(c(MAFA1,MAFA2,MAFB1,MAFB2),dim=c(1,4)))
        Lable     = rbind(Lable,paste(LA1,".",LA2,"+",LB1,".",LB2,sep=""))
      }
    }
    
  }
  
}
NComb = nrow(SNP_posi)

Power = array(0,dim=c(1,NComb))  #--- Power is used to calculate the power of the 5 methods
for (round in 1:NComb){
    cat(round,"\n")
    SNPs = SNP_posi[round,]

    for (i in 1:N.iter){
      
      #- simulate the data ------
      SData = RealisticSimulation(Sample_size=N,SNPs,Risk_Model=3, risk_effect=0.20)
      Y     = SData$Y
      X     = SData$X
      gene1 = SData$geno1
      gene2 = SData$geno2
      gene2 = gene2[,-14]
      gene2 = gene2[,-14]
      
      #-- Using Jung-Ying's Typical IBS -----
      Simi1     = JY.Typical.IBS(gene1,weight.power = -0.75)
      Simi2     = JY.Typical.IBS(gene2,weight.power = -0.75)
      
      pvalue    = JointTest(Y,X,Simi2,Simi1)
      
      if (pvalue<0.05) Power[1,round]=Power[1,round]+1/N.iter
    }
    
    #-- Output ----------
    #- The output is based on the SNP positions ---
}
final.output = cbind(SNP_posi,LDPatten, MAFPatten, Lable,t(Power))
write.table(final.output,"Model3.txt",sep=" ",row.name=F,col.name=F)
