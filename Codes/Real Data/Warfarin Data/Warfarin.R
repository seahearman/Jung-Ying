#----- This code is used to run the Warfarin Data

library(stats)
library(MASS)
library(CompQuadForm)
rm(list= ls())
set.seed(4)

#---------------------------- MAIN PART -------------------------------------

setwd(path.expand("D:\\Projects\\Jung-Ying\\Codes\\Main Functions"))

#---------------------- LOAD Functions -------------------------------------

source('TypicalIBS.r')
source('AverageIBS.r')
source('JointTest.r')
source('EpiTest.r')
source('OneGeneTest.r')

#---- Load data ----
setwd(path.expand("D:\\Projects\\Jung-Ying\\Codes\\Real Data\\Warfarin Data\\Raw data"))
Traits = as.matrix(read.table("Y.data",header=FALSE))
Y = array(Traits[,1], dim=c(nrow(Traits),1)) 
X = as.matrix(read.table("X.data",header=FALSE))

geno1 = as.matrix(read.table("genoA.data",header=FALSE))
Simi1 = AverageIBS_alleletype(geno1)
geno2 = as.matrix(read.table("genoB.data",header=FALSE))
Simi2 = TypicalIBS_alleletype(geno2)

Jointp = JointTest(Y,X,Simi1,Simi2)
Epip = EpiTest(Y,X,Simi1,Simi2)

Gene1p = OneGeneTest(Y,X,Simi2,Simi1)
Gene2p = OneGeneTest(Y,X,Simi1,Simi2)

cat("The Joint test for genoA and genoB is ",Jointp,"\n")
cat("The Epi test for genoA and genoB is ",Epip,"\n")
cat("The One gene test for gene1 is ", Gene1p, "\n")
cat("The One gene test for gene2 is ", Gene2p, "\n") 