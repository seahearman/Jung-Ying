# This program is used to convert the data into the pattern that my code can accept.
library(stats)
library(MASS)
rm(list=ls())
all.data <- read.table("gain.txt",header=FALSE)
dim_data  <- dim(all.data)
Nx3 <- dim_data[1]
N <- Nx3/3
new.data <- array(0,dim=c(N,4))

for (i in 1:4){
	all.temp <- array(all.data[,i],dim=c(Nx3,1))

	temp1 <- array(all.temp[seq(1,Nx3,by=3),1],dim=c(N,1))
	temp2 <- array(all.temp[seq(2,Nx3,by=3),1],dim=c(N,1))
	temp3 <-	array(all.temp[seq(3,Nx3,by=3),1],dim=c(N,1))
	
	temp <- temp1*0+ (temp2==1)+(temp3==1)*2
	new.data[,i]=temp

}

write.table(new.data,"gain2.txt",row.names=FALSE,col.names=FALSE)
 