library(stats)
library(MASS)
library(CompQuadForm)
library(pls)
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

AverageIBS_genotype = function (geno){
  #- If the data has a 2-allelic genotype format, i.e. there are only 2 allele for each
  #- marker, say A, a, and the genotype is recorded by a scale, that is, AA:0, Aa:1 and aa:2

  #- if the genotype has more than 2 allelic data, we may want to calculated in other way
  #- Sepetating the genotype into 2 haplotype matrix ----

  hA     		= geno				# hA is used to how many A's in the allele, it equals to genotype
	ha 	  		= 2-geno			# ha is used to how many a's in the allele, it equals to 2 - genotype
	
	N_gene  = ncol(geno)    # the number of markers
	S_temp  = hA %*% t(hA) + ha %*% t(ha)
	S         	=  S_temp/(4*N_gene)
  
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

#--------------------------------------------------------------------------------------
#                                       PCA                                           |
#--------------------------------------------------------------------------------------

 PCA_analysis <- function(phenotype,geno1,geno2){

  N   	<- nrow(phenotype)								# the nmuber of individuals in the sample
	
	#---------- 2.1 get the 1st PCA component of geno1 and geno2)----------------------

	gene1_PCA		 <- princomp(geno1,cor=FALSE, scores=TRUE)
	gene2_PCA		 <- princomp(geno2,cor=FALSE, scores=TRUE)
	
	Z1				 <- gene1_PCA$score[,1]
	Z2				 <- gene2_PCA$score[,1]
	
  dat1  	 <- data.frame(phenotype=phenotype,geno1,geno2,Z1*Z2)
	model_analysis 	 <- glm(phenotype~.,data = dat1)		# the alternative model
	
	ptest1		 	<- anova(model_analysis,test="Chisq")
	power1 		 	<- 1-pchisq((max(ptest1[[4]])-min(ptest1[[4]])),df=max(ptest1[[3]])-min(ptest1[[3]]))
	return(power1)
}
 
#--------------------------------------------------------------------------------------
#                                       PLS                                           |
#--------------------------------------------------------------------------------------
 
 PLS_analysis <- function(phenotype,geno1,geno2){

  N 		<- nrow(phenotype)								# the nmuber of individuals in the sample
	
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
	power1 		 <- 1-pchisq((max(ptest1[[4]])-min(ptest1[[4]])),df=max(ptest1[[3]])-min(ptest1[[3]]))
	
	return(power1)
}

LM_analysis <- function(Y,X,gene1,gene2){
  
    Inter = NULL
    
    for (i in 1:(ncol(gene1)-1)){
        for (j in i:(ncol(gene2))){
            Inter =  cbind(Inter,gene1[,i]*gene2[,j])
        }
    }   
    
    dat.all = data.frame(Y,X,gene1,gene2,Inter)
    fit0 = lm(Y~X,data = dat.all)     
    fit.all = update(fit0, . ~ . + gene1+gene2+Inter)
    p.all = anova(fit0,fit.all,test="F")
    return(p.all[2,6])
}
#---------------------------- Main part -------------------------------------

N     = 300
N.iter= 100

GeneA = array(c(1,3,7),dim=c(1,3))
LableA = array(c("HM","LC","LR"),dim=c(1,3))
GeneB = array(c(7,9,1,6,3),dim=c(1,5))
LableB = array(c("HC","NC","NM","LM","LR"),dim=c(1,5))

#--- Generate all the possible SNP combinations for Model 1 --
SNP_posi = NULL
Lable = NULL
for (gA1 in 1:2){
  
  SA1 = GeneA[1,gA1]
  LA1 = LableA[1,gA1]
  
  for (gA2 in (gA1+1):3){
    
    SA2 = GeneA[1,gA2]
    LA2 = LableA[1,gA2]
    
    for (gB1 in 1:4){
      
      SB1 = GeneB[1,gB1]
      LB1 = LableB[1,gB1]
      
      for (gB2 in (gB1+1):5){
        
        SB2 = GeneB[1,gB2]
        LB2 = LableB[1,gB2]
        
        SNP_posi = rbind(SNP_posi,array(c(SA1,SA2,SB1,SB2),dim=c(1,4)))
        Lable =rbind(Lable,paste(LA1,".",LB1,"+",LA2,".",LB2,sep=""))
        SNP_posi = rbind(SNP_posi,array(c(SA1,SA2,SB2,SB1),dim=c(1,4)))
        Lable =rbind(Lable,paste(LA1,".",LB2,"+",LA2,".",LB1,sep=""))
      }
    }
    
  }
  
}

NComb = nrow(SNP_posi)

Power = array(0,dim=c(5,NComb))  #--- Power is used to calculate the power of the 5 methods
for (round in 1:NComb){
    cat(round,"\n")
    SNPs = SNP_posi[round,]

    for (i in 1:N.iter){
      
      #- simulate the data ------
      SData = RealisticSimulation(Sample_size=N,SNPs,Risk_Model=1, risk_effect=0.15)
      Y     = SData$Y
      X     = SData$X
      gene1 = SData$geno1
      gene2 = SData$geno2
      
      #-- Using the Average IBS -----
      
      Simi1     = AverageIBS_genotype(gene1)
      Simi2     = AverageIBS_genotype(gene2)
      
      pvalue    = JointTest(Y,X,Simi2,Simi1)
      
      if (pvalue<0.05) Power[1,round]=Power[1,round]+1/N.iter
      
      #-- Using the Typical IBS -----
      
      Simi1     = TypicalIBS_genotype(gene1)
      Simi2     = TypicalIBS_genotype(gene2)
      
      pvalue    = JointTest(Y,X,Simi2,Simi1)
      
      if (pvalue<0.05) Power[2,round]=Power[2,round]+1/N.iter
      
          
      #-- Using the PCA -------
      
      pvalue = PCA_analysis(Y,gene1,gene2)
      if (pvalue<0.05) Power[3,round]=Power[3,round]+1/N.iter

      
      #-- Using the PLS -------
      pvalue = PLS_analysis(Y,gene1,gene2)
      if (pvalue<0.05) Power[4,round]=Power[4,round]+1/N.iter
      
      #-- Using the LM--------
      pvalue = LM_analysis(Y,X,gene1,gene2)
      if (pvalue<0.05) Power[5,round]=Power[5,round]+1/N.iter
    }
    
    #-- Output ----------
    #- The output is based on the SNP positions ---
}
final.output = cbind(Lable,t(Power))
write.table(final.output,"Model1.txt",sep=" ",row.name=F,col.name=F)

