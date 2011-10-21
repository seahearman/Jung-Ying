TypicalIBS_genotype  	 = function (geno){
  
  #- If the data has a 2-allelic genotype format, i.e. there are only 2 allele for each
  #- marker, say A, a, and the genotype is recorded by a scale, that is, AA:0, Aa:1 and aa:2

  #- if the genotype has more than 2 allelic data, we may want to calculated in other way
  #- Sepetating the genotype into 2 haplotype matrix ----
  #---------------------------------------------------------
  # this is the version can be applied to the 2-allelic data
  #---------------------------------------------------------
	N_gene = ncol(geno)
	N 		 = nrow(geno)

	h1A_ = array(0,dim=c(N,N_gene))			# h1A_ is used to indicate whether the 1st allele is A for gene
	h1_a = array(0,dim=c(N,N_gene))			# h1_a is used to indicate whether the 1st allele is a for gene
	h2A_ = array(0,dim=c(N,N_gene))			# h2A_ is used to indicate whether the 2nd allele is A for gene
	h2_a = array(0,dim=c(N,N_gene))			# h2_a is used to indicate whether the 2nd allele is a for gene

	h1_a = (geno>0)*1						# This is a small algorithm to generate the h1A_~h2_a from the genotype
	h2_a = geno-h1_a

	h1A_ = 1-h1_a
	h2A_ = 1-h2_a


	#----------------- Get the allele freq for each locus in each two gene-------------
	qA 	= array(1,dim=c(1,N_gene))		# calculate the freq of main  allele for gene
	qa	= array(1,dim=c(1,N_gene))		# calculate the freq of minor allele for gene



	wA 	=qA^(-0.75)			# using the freq^(-3/4) as the weight
	wa 	=qa^(-0.75)

	S_original = array(0,dim=c(N,N))						# S by {0, 0.5, 1} value set


	# Following is small loop to generate the similarity matrix for gene1 and gene2	

	for (i in 1:N_gene){
		temp1 = wA[i]*h1A_[,i] %*% t(h1A_[,i])+wa[i]*h1_a[,i] %*% t(h1_a[,i]) + wA[i]*h2A_[,i] %*% t(h2A_[,i])+wa[i]*h2_a[,i] %*% t(h2_a[,i])
		temp2 = wA[i]*h1A_[,i] %*% t(h2A_[,i])+wa[i]*h1_a[,i] %*% t(h2_a[,i]) + wA[i]*h2A_[,i] %*% t(h1A_[,i])+wa[i]*h2_a[,i] %*% t(h1_a[,i])
		S_original = S_original+ temp1*(temp1>=temp2)+temp2*(temp1<temp2)
	}

	S = S_original/(2*N_gene)
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

TypicalIBS_alleletype     = function (geno){
    
    #- if the marker are multi-allelic, then the genotype can not be recorded by just one value, it need to use a 2*1 vector to record
    #- the genotype. For this function, if there are M markers, 1*2M vector is needed for each individual.
    #--- how to generate the S matrix when the genotype is multi-allelic
    #---------------------------------------------------------
    # This is the version can be applied to multi-allelic data
    #---------------------------------------------------------
    N = nrow(geno)
    M = ncol(geno)
    S = array(0,dim=c(N,N))
  
    for (i in 1:N){
      for (j in 1:N){
        
        #- compute S[i,j]
        p1 = 1
        p2 = p1+1
        
        i1 = min(geno[i,p1],geno[i,p2])
        i2 = max(geno[i,p1],geno[i,p2])
        j1 = min(geno[j,p1],geno[j,p2])
        j2 = max(geno[j,p1],geno[j,p2])
        
        while (p2<=M){
          
          #---- if marker[i] is identical to marker[j]
          if ( (i1==j1) && (i2==j2))
            S[i,j] = S[i,j] + 2/M
          #---- if marker[i] totally differs from marker[j]
          else if ((i1!=j1) && (i1!= j2) && (i2!=j1) && (i2!=j2))
            S[i,j] = S[i,j]
          #---- if marker[i] share partially with marker[j] 
          else S[i,j] = S[i,j] + 1/M
         
          p1 = p2+1
          p2 = p1+1
        }
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