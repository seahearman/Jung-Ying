rm(list=ls())

n = 4
n.comb = n+choose(n,2)
A1 = array(0,dim=c(1,n.comb))
A2 = array(0,dim=c(1,n.comb))
p = 0
for (i in 1:n){
  for(j in i:n){
    p  = p+1
    A1[p] = i
    A2[p] = j
  }
}

S.average = array(0,dim=c(n.comb,n.comb))

for (i in 1:n.comb){
  for (j in 1:n.comb){
    if ((A1[i]==A2[i]) & (A1[j]==A2[j]) & (A1[i]==A1[j])){
      S.average[i,j] = 2
    }
    else if (A1[i]!= A1[j] & A1[i]!=A2[j] & A2[i]!=A1[j] & A2[i]!=A2[j] ){
      S.average[i,j] = 0
    }
  
  }
}

e1 = eigen(S.average)$values

S.typical = S.average

for (i in 1:n.comb){
  S.typical[i,i] = 2 
}

e2 = eigen(S.typical)$values
cat("average IBS: ",min(e1),"\n")
cat("typical IBS: ",min(e2),"\n")