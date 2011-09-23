# This program is used to convert the data into the pattern that my code can accept.
library(stats)
library(MASS)
rm(list=ls())

setwd(path.expand("D:/Written Prelim/Jung-Yng Zeng/data resort/Gene raw"))
AllGenes  	<-  list.files("D:/Written Prelim/Jung-Yng Zeng/data resort/Gene raw")
NoGene 	  	<-  length(AllGenes)

for (j in 1:NoGene){
		all.data 	<- as.matrix(read.table(AllGenes[j],header=FALSE))
		dim_data  	<- dim(all.data)
		Nx3 				<- dim_data[1]
		N 					<- Nx3/3
		N_gene 		<- dim_data[2] 
		new.data 	<- array(0,dim=c(N,N_gene))

		for (i in 1:N_gene){
		
				all.temp <- array(all.data[,i],dim=c(Nx3,1))

				temp1 <- array(all.temp[seq(1,Nx3,by=3),1],dim=c(N,1))
				temp2 <- array(all.temp[seq(2,Nx3,by=3),1],dim=c(N,1))
				temp3 <-	array(all.temp[seq(3,Nx3,by=3),1],dim=c(N,1))
				
				temp <- temp1*0+ (temp2==1)+(temp3==1)*2
				new.data[,i]=temp

		}
		output <- paste("D:/Written Prelim/Jung-Yng Zeng/data resort/Gene/",AllGenes[j],sep="")
		write.table(new.data, output, row.names=FALSE, col.names=FALSE)
}
 