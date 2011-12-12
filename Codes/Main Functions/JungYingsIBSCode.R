JY.Typical.IBS = function(SNP=addsnp, weight.power = (-1), Gmode="add", MAF=NULL)

{

  ## tSNP is SNPs-rows, individuals=columns format.
  
  tSNP = t(SNP)
  
  nsubj = ncol(tSNP)
  
  nloci = nrow(tSNP)
  
  if(Gmode=="add")
  
  {
  
    ## Weights for the S matrix calculated from the allele frequencies.
    
    tSNP.A = tSNP
    
    tSNP.a = 2 - tSNP.A
    
    AFA = rowSums(tSNP.A) / (2 * nsubj)
    
    AFa = 1 - AFA # Allele Freq of "a"
    
    wt.A = AFA ^ weight.power;
    
    wt.a = AFa ^ weight.power;
    
    ## use the following if use Wu weights
    
    ## wt.a = wt.A
    
    wttSNP.A = sqrt(wt.A) * tSNP.A
    
    wttSNP.a = sqrt(wt.a) * tSNP.a
    
    ## This is the avgIBS weighted S matrix.
    
    SS.wt = ( crossprod(wttSNP.A) + crossprod(wttSNP.a) ) / nloci
    
    ##-------------------
    
    ## SSmat by maxIBS kernel
    
    ##-------------------
    
    ## The maxIBS values can be calculated from the avgIBS S matrix
    
    ## plus an extra term for the het-het cases.
    
    mg2 = 1 - (1-tSNP)^2
    
    wtmg2.A = sqrt(wt.A) * mg2
    
    wtmg2.a = sqrt(wt.a) * mg2
    
    SS.hh = (crossprod(wtmg2.A) + crossprod(wtmg2.a)) / nloci
    
    SS.maxibs = SS.wt + SS.hh
    
    S      = SS.maxibs
    eg     =  eigen(S, symmetric=T)
    evalue =  eg$values
    le     =  evalue[1] # Biggest eigenvalue.
    if (le == 0){
      le   	= 1
  	}
  	tmpkey = (evalue / le) > 1e-7
  	l.maxibs = sum(tmpkey) # Rank of SSS
  
  }
  
  return(list("S"=SS.maxibs,"L"=l.maxibs))

}

JY.Average.IBS = function(SNP=addsnp, weight.power = (-1), Gmode="add", MAF=NULL)

{

  ## tSNP is SNPs-rows, individuals=columns format.
  
  tSNP = t(SNP)
  
  nsubj = ncol(tSNP)
  
  nloci = nrow(tSNP)
  
  if(Gmode=="add")
  
  {
  
    ## Weights for the S matrix calculated from the allele frequencies.
    
    tSNP.A = tSNP
    
    tSNP.a = 2 - tSNP.A
    
    AFA = rowSums(tSNP.A) / (2 * nsubj)
    
    AFa = 1 - AFA # Allele Freq of "a"
    
    wt.A = AFA ^ weight.power;
    
    wt.a = AFa ^ weight.power;
    
    ## use the following if use Wu weights
    
    ## wt.a = wt.A
    
    wttSNP.A = sqrt(wt.A) * tSNP.A
    
    wttSNP.a = sqrt(wt.a) * tSNP.a
    
    ## This is the avgIBS weighted S matrix.
    
    SS.wt = ( crossprod(wttSNP.A) + crossprod(wttSNP.a) ) / nloci
    
    S      = SS. wt
    eg     =  eigen(S, symmetric=T)
    evalue =  eg$values
    le     =  evalue[1] # Biggest eigenvalue.
    if (le == 0){
      le   	= 1
  	}
  	tmpkey = (evalue / le) > 1e-7
  	l.wt  	 		= sum(tmpkey) # Rank of SSS
  
  }
  
  return(list("S"=SS.wt,"L"=l.wt))

}