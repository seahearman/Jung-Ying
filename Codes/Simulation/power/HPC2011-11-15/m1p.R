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
for (gA1 in 1:2){
  
  SA1 = GeneA[1,gA1]
  LDA1 = ALD[1,gA1]
  MAFA1 = AMAF[1,gA1]
  LA1 = LableA[1,gA1]
  
  for (gA2 in (gA1+1):3){
    
    SA2 = GeneA[1,gA2]
    LA2 = LableA[1,gA2]
    LDA2 = ALD[1,gA2]
    MAFA2 = AMAF[1,gA2]
    
    for (gB1 in 1:4){
      
      SB1 = GeneB[1,gB1]
      LB1 = LableB[1,gB1]
      LDB1 = BLD[1,gB1]
      MAFB1 = BMAF[1,gB1]
      
      for (gB2 in (gB1+1):5){
        
        SB2 = GeneB[1,gB2]
        LB2 = LableB[1,gB2]
        LDB2 = BLD[1,gB2]
        MAFB2 = BMAF[1,gB2]
        
        SNP_posi = rbind(SNP_posi,array(c(SA1,SA2,SB1,SB2),dim=c(1,4)))
        LDPatten = rbind(LDPatten,array(c(LDA1,LDA2,LDB1,LDB2),dim=c(1,4)))
        MAFPatten = rbind(MAFPatten,array(c(MAFA1,MAFA2,MAFB1,MAFB2),dim=c(1,4)))
        Lable =rbind(Lable,paste(LA1,".",LB1,"+",LA2,".",LB2,sep=""))
        SNP_posi = rbind(SNP_posi,array(c(SA1,SA2,SB2,SB1),dim=c(1,4)))
        LDPatten = rbind(LDPatten,array(c(LDA1,LDA2,LDB2,LDB1),dim=c(1,4)))
        MAFPatten = rbind(MAFPatten,array(c(MAFA1,MAFA2,MAFB2,MAFB1),dim=c(1,4)))
        Lable =rbind(Lable,paste(LA1,".",LB2,"+",LA2,".",LB1,sep=""))
      }
    }
    
  }
  
}
final = cbind(SNP_posi,LDPatten, MAFPatten, Lable)
write.table(final,file="D:\\Projects\\Jung-Ying\\Codes\\Simulation\\power\\HPC2011-11-15\\model1 Patten.txt",sep=" ",row.name=F,col.name=F)
