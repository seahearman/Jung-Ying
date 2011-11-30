# This program is used for the EM algorithm.
library(stats)
library(MASS)
rm(list=ls())
set.seed(2)
#---------------------------- FUNCTIONS ----------------------------------------

#-------------------------FUNCTION 1: SIMULATION -------------------------------

  simulation = function(N, Tau, sigma){
    #1. Generate the fix effect
    X = array(1, dim=c(N,2))
    X[,2] = rnorm(N,mean=0, sd=1)
    
    Gamma = c(1,1)
    
    #2 Generate the random effect
    #2.1. Generate the H : haplotype matrix
    
    p = c(0.397,0.208,0.02,0.013,0.018,0.138,0.129,0.054)  # the propability of the 8 haplotypes in population
    H_all = t(rmultinom(N,2,p))
    
    #2.2. Generate beta    
    
    R_all = c( 1,0.8,0.8,0.6,0.6,0.2,0.8,0.6,
    0.8,1,0.6,0.8,0.4,0.4,0.6,0.8,
    0.8,0.6,1,0.8,0.8,0.4,0.6,0.4,
    0.6,0.8,0.8,1,0.6,0.6,0.4,0.6,
    0.6,0.4,0.8,0.6,1,0.6,0.4,0.2,
    0.2,0.4,0.4,0.6,0.6,1,0,0.2,
    0.8,0.6,0.6,0.4,0.4,0,1,0.8,
    0.6,0.8,0.4,0.6,0.2,0.2,0.8,1)   # the similarity matrix of the 8 haplotypes
    dim(R_all) = c(8,8)
 
    H=H_all
    R=R_all
    
    Sigma_beta=Tau*R
    
    beta=mvrnorm(n = 1, mu=rep(0,8), Sigma_beta, tol = 1e-6, empirical = FALSE)

    #3 Generate the error term
    
    e=rnorm(N,mean=0,sd=sqrt(sigma))
    
    #4. sum everything up!
    Y=X %*% Gamma+ H%*%beta+e
    return(list(Y = Y, X = X, H = H, R = R))
  }
  
  
#------------------------------- MAIN PART -------------------------------------  

  # basic data information
  
  N = 100 # The number of the objects
  True_Tau = 0.6
  True_sigma=1
  
  data_all=simulation(N,True_Tau,True_sigma) # Generate the simulated data
  
  Y=data_all$Y
  X=data_all$X
  H=data_all$H
  R=data_all$R
  
  p=ncol(X) # The number of fix effect
  q=ncol(H) # The number of haplotype
  P_X=X %*% solve(crossprod(X,X))%*% t(X) # The projection of X
  
  #----------------------- Initialize Tau and sigma ----------------------------
  
  Tau_new=0.6
  sigma_new=1
  
  #------------------------------- updata --------------------------------------
  Tau_record=array(0,dim=c(1,100))
  sigma_record=array(0,dim=c(1,100))
  for (i in 1:50){
    
    Tau=as.numeric(Tau_new)
#    Tau=True_Tau
    sigma=as.numeric(sigma_new)

    
    V = Tau*H%*%R%*%t(H)+sigma*diag(N)
#    inv_V=solve(V)
    inv_V=1/sigma*(diag(N)-Tau/sigma*H%*%R%*%solve(diag(q)+Tau/sigma*t(H)%*%H%*%R)%*%t(H))
    inv_VX=inv_V%*%X
    P = inv_V-inv_VX %*% solve(t(X)%*%inv_VX)%*%t(inv_VX)
    
    E_beta=Tau*R%*%t(H)%*%P%*%Y
    Var_beta=Tau*R-Tau^2*R%*%t(H)%*%P%*%H%*%R
    
    inv_R_E_beta=Tau*t(H)%*%P%*%Y
    inv_R_Var_beta=Tau*diag(q)-Tau^2*t(H)%*%P%*%H%*%R
    
    Tau_new=1/q*(t(E_beta)%*%inv_R_E_beta+sum(diag(inv_R_Var_beta)))
    sigma_new=1/(N-p)*(crossprod((Y-H%*%E_beta),(diag(N)-P_X))%*% (Y-H%*%E_beta)+sum(diag(crossprod(H, (diag(N)-P_X))%*%H%*%Var_beta)))
    Tau_record[i]=as.numeric(Tau_new)
    sigma_record[i]=as.numeric(sigma_new)
  }
#  Tau_record
#  sigma_record
  