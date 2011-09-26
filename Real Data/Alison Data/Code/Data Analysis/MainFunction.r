# Xin Wang, bioinformatics, NCSU, Aug, 2011 I use this to have a try

library(stats)
library(MASS)
library(CompQuadForm)
rm(list= ls())
set.seed(4)

#---------------------------- MAIN PART -------------------------------------  

setwd(path.expand("D:\\Projects\\Jung-Ying\\Real Data\\Alison Data\\Code\\Data Analysis"))

#----------------------  LOAD Functions -------------------------------------

source('TypicalIBS.r')
source('JointTest.r')
source('EpiTest.r')
source('OneGeneTest.r')

#---- Load data ----

Traits 		= as.matrix(read.table("Y.data",header=FALSE))
Y 				= array(Traits[,2], dim=c(301,1))		# Col 2 records BMI and Col 3 records log(BMI); it seems that using BMI would have a better result than log(BMI)
X					= as.matrix(read.table("X.data",header=FALSE))

geno1     = as.matrix(read.table("genoA.data",header=FALSE))
geno2     = as.matrix(read.table("genoB.data",header=FALSE))
Jointp    = JointTest(Y,X,geno1,geno2)
Epip      = EpiTest(Y,X,geno1,geno2)

Gene1p    = OneGeneTest(Y,X,geno2,geno1)
Gene2p    = OneGeneTest(Y,X,geno1,geno2)

cat("The Joint test for genoA and genoB is ",Jointp,"\n")
cat("The Epi test for genoA and genoB is ",Epip,"\n")
cat("The One gene test for gene1 is ", Gene1p, "\n")
cat("The One gene test for gene2 is ", Gene2p, "\n")  