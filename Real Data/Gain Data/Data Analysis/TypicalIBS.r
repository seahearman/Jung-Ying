#-------------- Function : generate the similarity matrix from the genotype-----------
# This funcion is used to calculate the typical IBS from the given genotype

TypicalIBS <- function (geno){

  #- Sepetating the genotype into 2 haplotype matrix ----

	hA 	  <- geno				# hA is used to how many A's in the allele, it equals to genotype
	ha 	  <- 2-geno			# ha is used to how many a's in the allele, it equals to 2 - genotype
	
	N_gene  <- ncol(geno)    # the number of markers
	S_temp  <- hA %*% t(hA) + ha %*% t(ha)
	S         <-  S_temp/(4*N_gene)
  
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
