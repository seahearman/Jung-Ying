TypicalIBS = function (geno){
  #- if the genotype has more than 2 allelic data, we may want to calculated in other way
  #- Sepetating the genotype into 2 haplotype matrix ----

#   hA 	  		= geno				# hA is used to how many A's in the allele, it equals to genotype
# 	ha 	  		= 2-geno			# ha is used to how many a's in the allele, it equals to 2 - genotype
# 	
# 	N_gene  = ncol(geno)    # the number of markers
# 	S_temp  = hA %*% t(hA) + ha %*% t(ha)
# 	S         	=  S_temp/(4*N_gene)
  
  #--- how to generate the S matrix when the genotype is multi-allelic
  
  N = nrow(geno)
  M = ncol(geno)
  
  S = array(0,dim=c(N,N))
  
  for (i in 1:N){
    for (j in 1:N){
      
      #- compute S[i,j]
      p1 = 1
      p2 = p1+1
      
      while (p2<=M){
        if ( (geno[i,p1]==geno[i,p2]) && (geno[i,p1]==geno[j,p1]) && (geno[j,p1]==geno[j,p2]))
          S[i,j] = S[i,j] + 2/M
        else if ((geno[i,p1]!=geno[j,p1]) && (geno[i,p1]!=geno[j,p2]) && (geno[i,p2]!=geno[j,p1]) && (geno[i,p2]!=geno[j,p2]))
          S[i,j] = S[i,j]
        else S[i,j] = S[i,j] + 1/M
        p1 = p2+1
        p2 = p1+1
      }
    }
      
  }
	
  eg     =  eigen(S, symmetric=T)
  evalue =  eg$values
  le     =  evalue[1] # Biggest eigenvalue.
	if (le == 0){
		le 		= 1
	}
	tmpkey = (evalue / le) > 1e-7
	ll  	 		= sum(tmpkey) # Rank of SSS
	RRRR	 	= diag(evalue[tmpkey])
	HHHH   	= eg$vectors[,tmpkey]
	return(list(S = S, R = RRRR, H = HHHH, L = ll))
	
}