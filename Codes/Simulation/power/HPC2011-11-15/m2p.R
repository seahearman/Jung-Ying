GeneA = array(c(1,3,7),dim=c(1,3))
LableA = array(c("HM","LC","LR"),dim=c(1,3))
GeneB = array(c(7,9,1,6,3),dim=c(1,5))
LableB = array(c("HC","NC","NM","LM","LR"),dim=c(1,5))
GeneA = array(c(1,3,7),dim=c(1,3))
ALD = array(c(0.588,0.126,0.159),dim=c(1,3))
AMAF = array(c(0.115,0.482,0.062),dim=c(1,3))
LableA = array(c("HM","LC","LR"),dim=c(1,3))
GeneB = array(c(7,9,1,6,3),dim=c(1,5))
BLD = array(c(0.259,0.23,0.186,0.136,0.065),dim=c(1,5))
BMAF = array(c(0.46,0.482,0.146,0.159,0.044),dim=c(1,5))
LableB = array(c("HC","NC","NM","LM","LR"),dim=c(1,5))

#--- Generate all the possible SNP combinations for Model 1 --
SNP_posi = NULL
Lable = NULL
LDPatten = NULL
MAFPatten = NULL

#--- Generate all the possible SNP combinations for Model 1 --
SNP_posi = NULL
Lable = NULL
for (gA1 in 1:3){
  for (gB1 in 1:5){
    
    SA1 = GeneA[1,gA1]
    LA1 = LableA[1,gA1]
    LDA1 = ALD[1,gA1]
    MAFA1 = AMAF[1,gA1]
    
    SB1 = GeneB[1,gB1]
    LB1 = LableB[1,gB1]
    LDB1 = BLD[1,gB1]
    MAFB1 = BMAF[1,gB1]
    
    SNP_posi = rbind(SNP_posi,array(c(SA1,SB1),dim=c(1,2)))
    Lable =rbind(Lable,paste(LA1,"+",LB1,sep=""))
    LDPatten = rbind(LDPatten,array(c(LDA1,LDB1),dim=c(1,2)))
    MAFPatten = rbind(MAFPatten,array(c(MAFA1,MAFB1),dim=c(1,2)))

  }  
}
final = cbind(SNP_posi,LDPatten, MAFPatten, Lable)
write.table(final,file="D:\\Projects\\Jung-Ying\\Codes\\Simulation\\power\\HPC2011-11-15\\model2 Patten.txt",sep=" ",row.name=F,col.name=F)
