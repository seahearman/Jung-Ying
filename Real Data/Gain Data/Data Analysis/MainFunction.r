# Xin Wang, bioinformatics, NCSU, Sep, 2011

library(stats)
library(MASS)
library(CompQuadForm)
rm(list= ls())
set.seed(4)

#--------------------------- Load functions ----------------------------------

setwd(path.expand("D:\\Projects\\Jung-Ying\\Real Data\\Gain Data\\Data Analysis"))
source("TypicalIBS.r")
source("JointTest.r")
source("EpiTest.r")
source("OneGeneTest.r")

#---------------------------- MAIN PART ------------------------------------- 

setwd(path.expand("D:\\Projects\\Jung-Ying\\Real Data\\Gain Data"))
	
#---- Load data ----

Traits   	<- as.matrix(read.table("good2_traits.txt",header=FALSE))
Y   			<- array(Traits[,2], dim=c(1068,1))		# Col 2 records BMI and Col 3 records log(BMI); it seems that using BMI would have a better result than log(BMI)
u 				<- array(1,dim=c(1068,1)) 						# Interuption
age 			<- array(Traits[,6], dim=c(1068,1))		# Col 6 records age
age2			<- age*age														# consider age^2 as a factor
sex				<- array(Traits[,5], dim=c(1068,1)) 	# Col 5 records sex
X					<- cbind(u,age,age2,sex)

setwd(path.expand("D:\\Projects\\Jung-Ying\\Real Data\\Gain Data\\Gene"))
AllGenes  <-  list.files("D:\\Projects\\Jung-Ying\\Real Data\\Gain Data\\Gene")
NoGene 	  <-  length(AllGenes)

cat("",file="D:\\Projects\\Jung-Ying\\Real Data\\Gain Data\\result.txt")

for (i in 2: (NoGene-1)){
	for (j in (i+1):NoGene){
	
		gene1 		<- as.matrix(read.table(AllGenes[i]))
		gene2 		<- as.matrix(read.table(AllGenes[j]))
    
    Combi   	<- paste(AllGenes[i],AllGenes[j])
		Combi		  <- gsub(".txt", " ",Combi)
    Combi  	  <- gsub(".TXT", " ",Combi)
    cat(Combi,"\n")
    
		Jointp 	  <- JointTest(Y,X,gene1,gene2)
    
    if (Jointp<0.1){
        Epip  <- EpiTest(Y,X,gene1,gene2)  
    }
    else{
        Epip  <- NA  
    }

    Gene1p    <- OneGeneTest(Y,X,gene2,gene1)
    Gene2p    <- OneGeneTest(Y,X,gene1,gene2)
		
    
		cat(Combi,Jointp,Epip,Gene1p,Gene2p,"\n",file="D:\\Projects\\Jung-Ying\\Real Data\\Gain Data\\result.txt",append=TRUE)
	}
}