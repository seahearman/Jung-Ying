setwd("D:\\Projects\\Jung-Ying\\Codes\\Simulation\\power\\HPC2011-11-15\\")
model= read.table("D:\\Projects\\Jung-Ying\\Codes\\Simulation\\power\\HPC2011-11-15\\Modle1.txt")
data=as.matrix(model)
m1=array(data,dim=c(60,5))
model= read.table("D:\\Projects\\Jung-Ying\\Codes\\Simulation\\power\\HPC2011-11-15\\Modle2.txt")
data=as.matrix(model)
m2=array(data,dim=c(15,5))
model= read.table("D:\\Projects\\Jung-Ying\\Codes\\Simulation\\power\\HPC2011-11-15\\Modle3.txt")
data=as.matrix(model)
m3=array(data,dim=c(30,5))
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

Lable1=Lable

SNP_posi = NULL
Lable = NULL
for (gA1 in 1:3){
  for (gB1 in 1:5){
    
    SA1 = GeneA[1,gA1]
    SA2 = 1
    LA1 = LableA[1,gA1]
    
    SB1 = GeneB[1,gB1]
    SB2 = 0
    LB1 = LableB[1,gB1]
    
    SNP_posi = rbind(SNP_posi,array(c(SA1,SA2,SB1,SB2),dim=c(1,4)))
    Lable =rbind(Lable,paste(LA1,"+",LB1,sep=""))

  }  
}
Lable2=Lable
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
        Lable =rbind(Lable,paste(LA1,".",LA2,"+",LB1,".",LB2,sep=""))
      }
    }
    
  }
  
}
Lable3=Lable
out1 = cbind(Lable1, m1)
out2 = cbind(Lable2, m2)
out3 = cbind(Lable3, m3)
write.table(out1,"m1.txt",sep=" ",row.name=F,col.name=F)
write.table(out2,"m2.txt",sep=" ",row.name=F,col.name=F)
write.table(out3,"m3.txt",sep=" ",row.name=F,col.name=F)
