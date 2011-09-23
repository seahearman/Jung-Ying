# Xin Wang, bioinformatics, NCSU, Aug, 2011

library(stats)
library(MASS)
library("pls")
rm(list= ls())
set.seed(4)

#--------------------------------------------------------------------------------------------------
#														FUNCTIONS																		|
#--------------------------------------------------------------------------------------------------
PCA_analysis <- function(phenotype,geno1,geno2){

  N 		<- nrow(phenotype)								# the nmuber of individuals in the sample
	
	#---------- 2.1 get the 1st PCA component of geno1 and geno2)----------------------
  
  geno1=geno1[,-5]
	gene1_PCA		 <- princomp(geno1,cor=FALSE, scores=TRUE)
	gene2_PCA		 <- princomp(geno2,cor=FALSE, scores=TRUE)
	
	Z1				 <- gene1_PCA$score[,1]
	Z2				 <- gene2_PCA$score[,1]
  
	In= Z1*Z2
  geno=cbind(geno1,geno2)

	model_analysis 	 <- glm(phenotype~geno1+geno2+In)		# the alternative model
	
	ptest1		 	<- anova(model_analysis,test="Chisq")
	po1 		 	<- 1-pchisq((max(ptest1[[4]])-min(ptest1[[4]])),df=max(ptest1[[3]])-min(ptest1[[3]]))
	return(po1)
}

#------------------- Function for PLS ----------------------------

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
	ptest1	 <- anova(fit3,test="Chisq")
	po1 		 <- 1-pchisq((max(ptest1[[4]])-min(ptest1[[4]])),df=max(ptest1[[3]])-min(ptest1[[3]]))
	
	return(po1)
}
#---------------------------- MAIN PART -------------------------------------  

setwd(path.expand("D:\\Projects\\Jung-Ying\\Real Data\\Alison Data"))

#---- Load data ----

Traits 		= as.matrix(read.table("Y.data",header=FALSE))
Y 				= array(Traits[,2], dim=c(301,1))		# Col 2 records BMI and Col 3 records log(BMI); it seems that using BMI would have a better result than log(BMI)
X					= as.matrix(read.table("X.data",header=FALSE))


fit0      = glm(Y~X)
Y         = residuals(fit0)

geno1     = as.matrix(read.table("genoA.data",header=FALSE))
geno2     = as.matrix(read.table("genoB.data",header=FALSE))



#--------------------------- Debug ------------------------------------

#   phenotype = Y
# 
# 
#   N   	<- nrow(phenotype)								# the nmuber of individuals in the sample
# 	
# 	#---------- 2.1 get the 1st PCA component of geno1 and geno2)----------------------
#   
#   geno1=geno1[,-5]
# 	gene1_PCA		 <- princomp(geno1,cor=FALSE, scores=TRUE)
# 	gene2_PCA		 <- princomp(geno2,cor=FALSE, scores=TRUE)
# 	
# 	Z1				 <- gene1_PCA$score[,1]
# 	Z2				 <- gene2_PCA$score[,1]
#   
# 	In= Z1*Z2
#   geno=cbind(geno1,geno2)
# 
# 	model_analysis 	 <- glm(phenotype~geno1+geno2+In)		# the alternative model
# 	
# 	anova(model_analysis,test="F")
# 
#   XX= cbind(geno,In)
# 
#   model_analysis 	 <- glm(phenotype~XX)		# the alternative model
# 	
# 	anova(model_analysis,test="F")

#   phenotype=Y
#   N=nrow(phenotype)
#   geno1=geno1[,-5]
#   center_geno1 <- apply(geno1, 2, function(x)(x-mean(x))/sd(x))
#   
#   center_geno2 <- apply(geno2, 2, function(x)(x-mean(x))/sd(x))
# 	
# 	#------------------------ The following 4 lines are used to remove NaN ------------
# 	bad			 <- sapply(center_geno1[1,], function(x) all(is.nan(x)))
# 	center_geno1 <- center_geno1[,!bad]
# 	bad			 <- sapply(center_geno2[1,], function(x) all(is.nan(x)))
# 	center_geno2 <- center_geno2[,!bad]	
#   
# 
# 
# 	Y 			 <- (phenotype-mean(phenotype))/sd(phenotype)	
# 	
# 	dat1		 <- data.frame(phenotype=phenotype,center_geno2)
# 	fit1		 <- glm(phenotype~.,data=dat1)
# 	mu			 <- fitted.values(fit1)
# 	
# 	dat2 		 <- data.frame(Gene1=I(center_geno1),Y_Gene2=I(cbind(center_geno2,Y)))
# 	pls1		 <- plsr(Y_Gene2 ~ Gene1, data = dat2)
# 	In			 <- scores(pls1)[,1]*mu
# 	dat3		 <- data.frame(phenotype=phenotype,geno1,geno2,In)
# 	fit3		 <- glm(phenotype~.,data=dat3)
#   anova(fit3,test="F")
#   XX= cbind(geno1,geno2,In)
# 
#   fit4		 <- glm(phenotype~XX)
#   anova(fit4,test="F")





#-------------------------------------------------------------------------





# geno1     = as.matrix(read.table("genoA.data",header=FALSE))
# geno2     = as.matrix(read.table("genoB.data",header=FALSE))
# PCAp      = PCA_analysis(Y,geno1,geno2)
# PLSp      = PLS_analysis(Y,geno1,geno2)
# 
# cat("The Joint test for PCA is ",PCAp,"\n")
# cat("The Joint test for PLS is ",PLSp,"\n")