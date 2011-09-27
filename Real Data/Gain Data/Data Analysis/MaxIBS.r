MaxIBS  	 <- function (geno){

	N_gene <- ncol(geno)
	
	h1A_ <- array(0,dim=c(N,N_gene))			# h1A_ is used to indicate whether the 1st allele is A for gene
	h1_a <- array(0,dim=c(N,N_gene))			# h1_a is used to indicate whether the 1st allele is a for gene
	h2A_ <- array(0,dim=c(N,N_gene))			# h2A_ is used to indicate whether the 2nd allele is A for gene
	h2_a <- array(0,dim=c(N,N_gene))			# h2_a is used to indicate whether the 2nd allele is a for gene

	h1_a <- (geno>0)*1						# This is a small algorithm to generate the h1A_~h2_a from the genotype
	h2_a <- geno-h1_a

	h1A_ <- 1-h1_a
	h2A_ <- 1-h2_a


	#----------------- Get the allele freq for each locus in each two gene-------------
	qA 	<- array(1,dim=c(1,N_gene))		# calculate the freq of main  allele for gene
	qa		<- array(1,dim=c(1,N_gene))		# calculate the freq of minor allele for gene



	wA 	<-qA^(-0.75)			# using the freq^(-3/4) as the weight
	wa 	<-qa^(-0.75)

	S_original <- array(0,dim=c(N,N))						# S by {0, 0.5, 1} value set


	# Following is small loop to generate the similarity matrix for gene1 and gene2	

	for (i in 1:N_gene){
		temp1 <- wA[i]*h1A_[,i] %*% t(h1A_[,i])+wa[i]*h1_a[,i] %*% t(h1_a[,i]) + wA[i]*h2A_[,i] %*% t(h2A_[,i])+wa[i]*h2_a[,i] %*% t(h2_a[,i])
		temp2 <- wA[i]*h1A_[,i] %*% t(h2A_[,i])+wa[i]*h1_a[,i] %*% t(h2_a[,i]) + wA[i]*h2A_[,i] %*% t(h1A_[,i])+wa[i]*h2_a[,i] %*% t(h1_a[,i])
		S_original <- S_original+ temp1*(temp1>=temp2)+temp2*(temp1<temp2)
	}

	S <- S_original/(2*N_gene)
	
  for (i in 1:nrow(S)){
    S[i,i]=1
	}
  eg        <-  eigen(S, symmetric=T)
  evalue  <-  eg$values
  le        <-  evalue[1] # Biggest eigenvalue.
	if (le == 0){
		le <- 1
	}
	tmpkey <- (evalue / le) > 1e-7
	ll  	 <- sum(tmpkey) # Rank of SSS
	RRRR	 <- diag(evalue[tmpkey])
	HHHH   <- eg$vectors[,tmpkey]
	return(list(S = S, R = RRRR, H = HHHH, L = ll))	

}