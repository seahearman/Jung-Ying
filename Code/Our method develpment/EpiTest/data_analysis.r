# This program is used to run the real data.
library(stats)
library(MASS)
rm(list=ls())
set.seed(3)
# hsreg_VC <- function(Y,genotype,N_gene1,N_gene2){

phenotype <- read.table("good2_traits.txt",header=FALSE)
Y <- array(phenotype[,4],dim=c(1068,1))
genotype <- read.table("gain2.txt",header=FALSE)
N_gene1 <- 2
N_gene2 <- 2

	
N 		<- nrow(Y)						# the number of individual
X 		<- array(1,dim=c(N,1))			# the population mean

#------------------ Using genotype to get the S -----------------------------------

geno1 	<- genotype[,1:N_gene1]							# genotype for 1st gene
geno2 	<- genotype[,(N_gene1+1):(N_gene1+N_gene2)]		# genotype for 2nd gene

#- Sepetating the genotype into 2 haplotype matrix ----
#-- for the geno1 --
g1h1A_ <- array(0,dim=c(N,N_gene1))			# g1h1A_ is used to indicate whether the 1st allele is A for gene1
g1h1_a <- array(0,dim=c(N,N_gene1))			# g1h1_a is used to indicate whether the 1st allele is a for gene1
g1h2A_ <- array(0,dim=c(N,N_gene1))			# g1h2A_ is used to indicate whether the 2nd allele is A for gene1
g1h2_a <- array(0,dim=c(N,N_gene1))			# g1h2_a is used to indicate whether the 2nd allele is a for gene1

g1h1_a <- (geno1>0)*1						# This is a small algorithm to generate the g1h1A_~g1h2_a from the genotype
g1h2_a <- geno1-g1h1_a

g1h1A_ <- 1-g1h1_a
g1h2A_ <- 1-g1h2_a

#-- for the geno2 --
g2h1A_ <- array(0,dim=c(N,N_gene2))			# This is similar to gene1
g2h1_a <- array(0,dim=c(N,N_gene2))
g2h2A_ <- array(0,dim=c(N,N_gene2))
g2h2_a <- array(0,dim=c(N,N_gene2))

g2h1_a <- (geno2>0)*1
g2h2_a <- geno2-g2h1_a

g2h1A_ <- 1-g2h1_a
g2h2A_ <- 1-g2h2_a

#----------------- Get the allele freq for each locus in each two gene-------------
q1A 	<- array(1,dim=c(1,N_gene1))		# calculate the freq of main  allele for gene1
q1a		<- array(1,dim=c(1,N_gene1))		# calculate the freq of minor allele for gene1
q2A 	<- array(1,dim=c(1,N_gene2))		# calculate the freq of main  allele for gene2
q2a		<- array(1,dim=c(1,N_gene2))		# calculate the freq of minor allele for gene2


w1A 	<-q1A^(-0.75)			# using the freq^(-3/4) as the weight
w1a 	<-q1a^(-0.75)
w2A 	<-q2A^(-0.75)
w2a 	<-q2a^(-0.75)	

S1_original <- array(0,dim=c(N,N))						# S by {0, 0.5, 1} value set
S2_original <- array(0,dim=c(N,N))						# S by {0, 0.5, 1} value set

# Following is small loop to generate the similarity matrix for gene1 and gene2	

for (i in 1:N_gene1){
	temp1 <- w1A[i]*g1h1A_[,i] %*% t(g1h1A_[,i])+w1a[i]*g1h1_a[,i] %*% t(g1h1_a[,i]) + w1A[i]*g1h2A_[,i] %*% t(g1h2A_[,i])+w1a[i]*g1h2_a[,i] %*% t(g1h2_a[,i])
	temp2 <- w1A[i]*g1h1A_[,i] %*% t(g1h2A_[,i])+w1a[i]*g1h1_a[,i] %*% t(g1h2_a[,i]) + w1A[i]*g1h2A_[,i] %*% t(g1h1A_[,i])+w1a[i]*g1h2_a[,i] %*% t(g1h1_a[,i])
	S1_original <- S1_original+ temp1*(temp1>=temp2)+temp2*(temp1<temp2)
}

for (i in 1:N_gene2){
	temp1 <- w2A[i]*g2h1A_[,i] %*% t(g2h1A_[,i])+w2a[i]*g2h1_a[,i] %*% t(g2h1_a[,i]) + w2A[i]*g2h2A_[,i] %*% t(g2h2A_[,i])+w2a[i]*g2h2_a[,i] %*% t(g2h2_a[,i])
	temp2 <- w2A[i]*g2h1A_[,i] %*% t(g2h2A_[,i])+w2a[i]*g2h1_a[,i] %*% t(g2h2_a[,i]) + w2A[i]*g2h2A_[,i] %*% t(g2h1A_[,i])+w2a[i]*g2h2_a[,i] %*% t(g2h1_a[,i])
	S2_original <- S2_original+ temp1*(temp1>=temp2)+temp2*(temp1<temp2)
}

S1 <- S1_original/(2*N_gene1)			
S2 <- S2_original/(2*N_gene2)
S12 <- S1* S2

Q 		<- diag(N) - X %*% solve(t(X) %*% X) %*% t(X)

sigma 	<- as.numeric (t(Y) %*% Q %*% Y/(N-ncol(X)))

P0 		<- 1 / sigma * Q

T 		<- 1 / (2*sigma^2) * t(Y) %*% Q %*% (S1 + S2 + S12) %*% Q %*% Y
Weights_all <- eigen(1 / (2 * sigma) * Q %*% (S1 + S2 + S12) %*% Q, symmetric=TRUE, only.value=TRUE)
temp 	<- Weights_all$values
temp2 	<- sort(temp,decreasing=TRUE)
dim(temp2) <- c(N,1)
big_enough <- sum(temp>10^-3)			# Get the number of big eigen values. here, the threshold for "big" is 10^-3 

Weights <- array(temp2[1:big_enough,1],dim=c(big_enough,1))


C 		<- rchisq(n=big_enough*5000,df=1)
dim(C) 	<- c(5000,big_enough)	# so now the dimension of C is L*5000

weightC <- C %*% Weights 
N_w 	<- nrow(weightC)
p 		<- sum(as.numeric(T)<=weightC)/N_w
	


