#--------------------------------------------------------------------------------------
#                                       PCA                                           |
#--------------------------------------------------------------------------------------

 PCA_analysis <- function(phenotype,geno1,geno2){

  N 		<- nrow(phenotype)								# the nmuber of individuals in the sample
	
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
