
 # This program is used for the EM algorithm.
library(stats)
library(MASS)
rm(list= ls())
set.seed(2)
 #---------------------------- FUNCTIONS ----------------------------------------

 #-------------------------FUNCTION 1: SIMULATION -------------------------------

 # This function is to genearate data accroding to the model Y=X*gammar+H^A * beta^A + H^B * beta^B + e

simulation <- function(N, Tau_A,Tau_B,phi, sigma){
	#1. Generate the fix effect
	X <- array(1, dim <- c(N,2))
	X[,2] <- rnorm(N,mean=0, sd=1)

	Gamma <- c(1,1)

	#2 Generate the random effect
	#2.1. Generate the H_A and H_B : haplotype matrix

	p <- c(0.397,0.208,0.02,0.013,0.018,0.138,0.129,0.054)  # the propability of the 8 haplotypes in population
	H_all <- t(rmultinom(N,2,p))
	H_A <- H_all

	H_all <- t(rmultinom(N,2,p))
	H_B <- H_all

	
	qA <- length(p)
	qB <- length(p)
	
	H_AB <- array(0, dim <- c(N,qA*qB))
	# 2.2 Genearte the H_AB from H_A and H_B
	pcol <- 0
	for (A in 1:qA){
		for (B in 1:qB){
			pcol <- pcol+1
			H_AB[, pcol] <- H_A[,A]*H_B[,B]
		}
	}

	#2.3 Generate the variance of beta: R_A and R_B    

	R_all <- c( 1,0.8,0.8,0.6,0.6,0.2,0.8,0.6,
			 0.8,1,0.6,0.8,0.4,0.4,0.6,0.8,
			 0.8,0.6,1,0.8,0.8,0.4,0.6,0.4,
			 0.6,0.8,0.8,1,0.6,0.6,0.4,0.6,
			 0.6,0.4,0.8,0.6,1,0.6,0.4,0.2,
			 0.2,0.4,0.4,0.6,0.6,1,0,0.2,
			 0.8,0.6,0.6,0.4,0.4,0,1,0.8,
			 0.6,0.8,0.4,0.6,0.2,0.2,0.8,1)   # the similarity matrix of the 8 haplotypes

	dim(R_all) <- c(8,8)

	R_A <- R_all
	R_B <- R_all

	
	#2.4 generate R_AB according to R_A and R_B (there are 4 loops, really ineffiecent)
	R_AB <- array(0, dim <- c(qA*qB,qA*qB))
	prow <- 0
	for (A1 in 1:qA){
		for (B1 in 1:qB){
			prow <- prow+1
			pcol <- 0
			for (A2 in 1:qA){
				for (B2 in 1:qB){
					pcol <- pcol+1
					R_AB[prow, pcol] <- R_A[A1,A2]*R_B[B1,B2]
				}
			}
		}
	}
	Sigma_beta_A <- Tau_A*R_A
	Sigma_beta_B <- Tau_B*R_B
	Sigma_beta_AB <- phi*R_AB

	beta_A <- mvrnorm(n = 1, mu=rep(0,qA), Sigma_beta_A, tol = 1e-6, empirical = FALSE)
	beta_B <- mvrnorm(n = 1, mu=rep(0,qB), Sigma_beta_B, tol = 1e-6, empirical = FALSE)
	beta_AB <- mvrnorm(n = 1, mu=rep(0,qA*qB), Sigma_beta_AB, tol = 1e-6, empirical = FALSE)
	#3 Generate the error term

	e <- rnorm(N,mean=0,sd=sqrt(sigma))

	#4. sum everything up!
	Y <- X %*% Gamma+ H_A %*% beta_A + H_B%*%beta_B+H_AB %*% beta_AB+e
	return(list(Y = Y, X = X, H_A = H_A,H_B = H_B, H_AB=H_AB, R_A = R_A,R_B = R_B, R_AB=R_AB))
}


#---------------------------- MAIN PART -------------------------------------  

# basic data information

N <- 50 # The number of the objects
True_Tau_A <- 0.6
True_Tau_B <- 0.4
True_phi <- 0
True_sigma <- 1
TA <- array(0,dim <- c(1,100))
TB <- array(0,dim <- c(1,100))
S <- array(0,dim <- c(1,100))
a_record <- array(0,dim <- c(1,100))
b_record <- array(0,dim <- c(1,100))
T0_record <- array(0,dim <- c(1,100))
ii=1
for (ii in 1:100){
	data_all <- simulation(N,True_Tau_A,True_Tau_B,True_phi,True_sigma) # Generate the simulated data

	Y <- data_all$Y
	X <- data_all$X
	H_A <- data_all$H_A
	H_B <- data_all$H_B
	H_AB <- data_all$H_AB
	R_A <- data_all$R_A
	R_B <- data_all$R_B
	R_AB <- data_all$R_AB

	p <- ncol(X) # The number of fix effect
	q_A <- ncol(H_A) # The number of A haplotype
	q_B <- ncol(H_B) # The number of B haplotype
	P_X <- X %*% solve(crossprod(X,X))%*% t(X) # The projection of X

	#----------------------- Initialize Tau and sigma ----------------------------

	Tau_A_new <- 0.5                        
	Tau_B_new <- 0.5
	sigma_new <- 1

	#------------------------------- updata --------------------------------------
	Tau_record <- array(0,dim <- c(1,100))  
	sigma_record <- array(0,dim <- c(1,100))
	R <- array(0,dim <- c(q_A+q_B,q_A+q_B))
	for (i in 1:50){
		
		Tau_A<- as.numeric(Tau_A_new)
		Tau_B<- as.numeric(Tau_B_new)
		R <- rbind(cbind(Tau_A*R_A,array(0,dim <- c(q_A,q_B))),cbind(array(0,dim <- c(q_B,q_A)),Tau_B*R_B))
		H <- cbind(H_A,H_B)

		sigma <- as.numeric(sigma_new)


		V <- H%*%R%*%t(H)+sigma*diag(N)

		inv_V <- 1/sigma*(diag(N)-1/sigma*H %*% R %*% solve(diag(q_A+q_B)+1/sigma*t(H)%*%H%*%R)%*%t(H))
		inv_VX <- inv_V%*%X                   #    To speed up the computation
		P <- inv_V-inv_VX %*% solve(t(X)%*%inv_VX)%*%t(inv_VX)

		E_beta <- R %*% t(H) %*% P %*% Y
		E_betaA <- E_beta[1:q_A]
		E_betaB <- E_beta[(q_A+1):(q_A+q_B)]

		Var_beta <- R-R%*%t(H)%*%P%*%H%*%R
		Var_betaA <- Var_beta[1:q_A,1:q_A]
		Var_betaB <- Var_beta[(q_A+1):(q_A+q_B),(q_A+1):(q_A+q_B)] 

		inv_RE_betaA <- Tau_A*t(H_A)%*%P%*%Y
		inv_RVar_betaA <- Tau_A*diag(q_A)-Tau_A^2*t(H_A)%*%P%*%H_A%*%R_A

		inv_RE_betaB <- Tau_B*t(H_B)%*%P%*%Y
		inv_RVar_betaB <- Tau_B*diag(q_B)-Tau_B^2*t(H_B)%*%P%*%H_B%*%R_B

		Tau_A_new <- 1/q_A*(t(E_betaA)%*%inv_RE_betaA+sum(diag(inv_RVar_betaA)))
		Tau_B_new <- 1/q_B*(t(E_betaB)%*%inv_RE_betaB+sum(diag(inv_RVar_betaB)))
		sigma_new <- 1/(N-p)*(crossprod((Y-H%*%E_beta),(diag(N)-P_X))%*% (Y-H%*%E_beta)+sum(diag(crossprod(H, (diag(N)-P_X))%*%H%*%Var_beta)))
	}
	TA[ii] <- as.numeric(Tau_A_new)     #  Tau_record
	TB[ii] <- as.numeric(Tau_B_new)     #  Tau_record
	S[ii] <- as.numeric(sigma_new)

	# Get the new P0
	P0 <- array(0,dim <- c(N,N))
	P0 <- P                               #Here under the null hypothesis, P0=P

	S_A <- H_A %*% R_A %*% t(H_A)
	S_B <- H_B %*% R_B %*% t(H_B)
	S_AB <- H_AB %*% R_AB %*% t(H_AB)
	
	T0 <- 1/2*t(Y) %*% P0 %*% S_AB %*% P0 %*% Y
	E_T0 <- 1/2*sum(diag(P0 %*% S_AB))	#the expectation of T0
	I_pp <- 1/2*sum(diag(P0 %*% S_AB %*% P0 %*% S_AB))
	I_tp <- c(1/2*sum(diag(P0 %*% S_AB %*% P0 %*% S_A)), 1/2*sum(diag(P0 %*% S_AB %*% P0 %*% S_B)), 1/2*sum(diag(P0 %*% S_AB %*% P0 ))) 
	I_pt <- t(I_tp)
	I_tt <- array(0,dim <- c(3,3))
	I_tt[1,] <- cbind(1/2*sum(diag(P0 %*% S_A %*% P0 %*% S_A)), 1/2*sum(diag(P0 %*% S_A %*% P0 %*% S_B)), 1/2*sum(diag(P0 %*% S_A %*% P0 )))
	I_tt[2,] <- cbind(1/2*sum(diag(P0 %*% S_B %*% P0 %*% S_A)), 1/2*sum(diag(P0 %*% S_B %*% P0 %*% S_B)), 1/2*sum(diag(P0 %*% S_B %*% P0 )))
	I_tt[3,] <- cbind(1/2*sum(diag(P0 %*% P0 %*% S_A)), 1/2*sum(diag(P0 %*% P0 %*% S_B)), 1/2*sum(diag(P0 %*% P0 )))
	Var_T0 <- I_pp-I_pt %*% solve(I_tt) %*% I_tp
	b_record[ii] <- Var_T0/E_T0
	a_record[ii] <- E_T0/b_record[ii]
	T0_record[ii] <- T0/b_record[ii]
}                      
a_mean <- mean(a_record)
require(graphics)    
qqplot(T0_record, rgamma(200,a_mean))            