ModifiedIBS = function (geno){

  #- Sepetating the genotype into 2 haplotype matrix ----

	hA 	 		 = geno				# hA is used to how many A's in the allele, it equals to genotype
	ha 	  		= 2-geno			# ha is used to how many a's in the allele, it equals to 2 - genotype
	
	N 		  	= nrow(geno)    # sample size
	N_gene = ncol(geno)    # the number of markers
	
	S 		  	= array(0,dim=c(N,N))
	
	#--- Generate the Sij -----
	
	for (i in 1:N){
		for(j in i:N){
			
			tempS 				= (hA[i,]*hA[j,]+ha[i,]*ha[j,])/(4*N_gene)
			dim(tempS)  	=  c(1,N_gene)
			S[i,j] 				= sum(tempS)
			
			for (i1 in 1:(N_gene-1)){
				for (j1 in (i1+1):N_gene){
				
					S[i,j]		= S[i,j]+tempS[1,i1]*tempS[1,j1]
					
				}
			}
			
			S[j,i] 				= S[i,j]
		}
	}
    eg        =  eigen(S, symmetric=T)
    evalue  =  eg$values
    le        =  evalue[1] # Biggest eigenvalue.
	if (le == 0){
		le = 1
	}
	tmpkey = (evalue / le) > 1e-7
	ll  	 = sum(tmpkey) # Rank of SSS
	RRRR	 = diag(evalue[tmpkey])
	HHHH   = eg$vectors[,tmpkey]
	return(list(S = S, R = RRRR, H = HHHH, L = ll))
	
}