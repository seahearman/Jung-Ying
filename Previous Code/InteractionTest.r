# Xin Wang, bioinformatics, NCSU, Oct, 2010

 #This is a simulation program to sovle the 2 Genes (including the interaction) modle.
 #Here I first simulate the data under the null hypothesis: phi=0, then use the score test to see how is the statsitcs behave.
 #When we only have one gene in the model, the approximated gamma distribution can describe the Ti statistics very well, 
 #but since now we have a more complicated model, we need to estimate tauA and tauB (sigma is always estimated whatever model),
 #it turns out the estimates of the two taus are not very good (especially when tau is small) which affect the approximated gamma
 #distribution. So finally we decide to try the sum of weighted chi-square distriubion to see how this distribution works.  

library(stats)
library(MASS)
rm(list= ls())
set.seed(4)
 #---------------------------- FUNCTIONS ----------------------------------------

 #-------------------------FUNCTION 1: SIMULATION -------------------------------

 # The simulation function is to genearate data accroding to the model Y=X*gammar+H^A * beta^A + H^B * beta^B + e

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

#--------------- Function 2: weighted chi-square distribution----------------
# The main idea for for this function is:
# 1. generate 8*10,000 chi square distributed random variables: Matrix C(this should be outside of the function to save time,
# you don't need to generate lots of random vairables each time)
# 2. Acoording to the eigenvalues, calculated 10,000 weighted chi-square distributed random variables
# 3. given the Ti statistics, get the p-value.

weighted_chi2 <- function(eigenvalues,Ti){
	c1 <- sum(eigenvalues)
	c2 <- sum(eigenvalues^2)
	c3 <- sum(eigenvalues^3)
	h <- c2^3/c3^2
	y <- (Ti-c1)*(h/c2)^(0.5)+h
	p <- pchisq(y,df=h)
	return(p_value=p)
}

#---------------------------- MAIN PART -------------------------------------  

# basic data information
Run_parameter <- function(N,True_Tau_A,True_Tau_B){
	True_phi   	<- 0									# True value for phi which under the null hypothesis, phi usually equals 0
	True_sigma 	<- 1									# True value for sigma
	round_all  	<- 500  								# The total times we will run the simulations
	T0_record 	<- array(0,dim <- c(1,round_all))		# The array used to record the T0 statistics for each simulation
	P_chi		<- array(0,dim <- c(1,round_all)) 		# The array used to record the p_value from the weighted chi-square distribution for each simulation
	tau_A_record <- array(0,dim <- c(1,round_all))
	tau_B_record <- array(0,dim <- c(1,round_all))

	for (ii in 1:round_all){
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

		Tau_A_new <- True_Tau_A                        
		Tau_B_new <- True_Tau_B
		sigma_new <- True_sigma

		#------------------------------- updata --------------------------------------

		R <- array(0,dim <- c(q_A+q_B,q_A+q_B))
		repeat{		
			Tau_A	<- as.numeric(Tau_A_new)
			Tau_B	<- as.numeric(Tau_B_new)
			R 			<- rbind(cbind(Tau_A*R_A,array(0,dim <- c(q_A,q_B))),cbind(array(0,dim <- c(q_B,q_A)),Tau_B*R_B))
			H 			<- cbind(H_A,H_B)

			sigma 	<- as.numeric(sigma_new)


			V 								<- H%*%R%*%t(H)+sigma*diag(N)

			inv_V 						<- 1/sigma*(diag(N)-1/sigma*H %*% R %*% solve(diag(q_A+q_B)+1/sigma*t(H)%*%H%*%R)%*%t(H))
			inv_VX 					<- inv_V%*%X                   #    To speed up the computation
			P 								<- inv_V-inv_VX %*% solve(t(X)%*%inv_VX)%*%t(inv_VX)

			E_beta 					<- R %*% t(H) %*% P %*% Y
			E_betaA 					<- E_beta[1:q_A]
			E_betaB 					<- E_beta[(q_A+1):(q_A+q_B)]

			Var_beta				<- R-R%*%t(H)%*%P%*%H%*%R
			Var_betaA 			<- Var_beta[1:q_A,1:q_A]
			Var_betaB 			<- Var_beta[(q_A+1):(q_A+q_B),(q_A+1):(q_A+q_B)] 

			inv_RE_betaA 		<- Tau_A*t(H_A)%*%P%*%Y
			inv_RVar_betaA <- Tau_A*diag(q_A)-Tau_A^2*t(H_A)%*%P%*%H_A%*%R_A

			inv_RE_betaB		<- Tau_B*t(H_B)%*%P%*%Y
			inv_RVar_betaB <- Tau_B*diag(q_B)-Tau_B^2*t(H_B)%*%P%*%H_B%*%R_B

			Tau_A_new 			<- 1/q_A*(t(E_betaA)%*%inv_RE_betaA+sum(diag(inv_RVar_betaA)))
			Tau_B_new 			<- 1/q_B*(t(E_betaB)%*%inv_RE_betaB+sum(diag(inv_RVar_betaB)))
			sigma_new 			<- 1/(N-p)*(crossprod((Y-H%*%E_beta),(diag(N)-P_X))%*% (Y-H%*%E_beta)+sum(diag(crossprod(H, (diag(N)-P_X))%*%H%*%Var_beta)))
			
			diff_A 					<- abs((as.numeric(Tau_A_new)-Tau_A)/as.numeric(Tau_A_new))		# use the relative difference rather than the absoult difference
			diff_B 					<- abs((as.numeric(Tau_B_new)-Tau_B)/as.numeric(Tau_B_new))
			diff_s 					<- abs((as.numeric(sigma_new)-sigma)/as.numeric(sigma_new))
			if ((diff_A<0.001) & (diff_B<0.001) & (diff_s<0.001)) break
		}
		tau_A_record[1,ii] <- as.numeric(Tau_A_new)-True_Tau_A
		tau_B_record[1,ii] <- as.numeric(Tau_B_new)-True_Tau_B
		# Get the new P0
		P0 <- array(0,dim <- c(N,N))
		P0 <- P                               			#Here under the null hypothesis, P0=P
		
		S_A <- H_A %*% R_A %*% t(H_A)
		S_B <- H_B %*% R_B %*% t(H_B)
		S_AB <- H_AB %*% R_AB %*% t(H_AB)
		
		T0 <- 1/2*t(Y) %*% P0 %*% S_AB %*% P0 %*% Y	#Get T0 under null hypothesis
		e <- eigen(V, symmetric=TRUE)
		V_eigen <- e$vectors
		V_square_root <- V_eigen %*% diag(sqrt(e$values)) %*% t(V_eigen)


		Weights_all <- eigen(1/2*V_square_root %*% P0 %*% S_AB %*% P0 %*% V_square_root, symmetric=TRUE, only.value=TRUE)
		temp <- Weights_all$values
		temp2 <- sort(temp,decreasing=TRUE)
		dim(temp2) <- c(N,1)
		Weights <- array(temp2[1:(q_A*q_B),1],dim=c(q_A*q_B,1))
		P_chi[ii] <- weighted_chi2(Weights,T0)
	} 
	biasA <- mean(tau_A_record)
	varA <- mean(tau_A_record^2)
	biasB <- mean(tau_B_record)
	varB <- mean(tau_B_record^2)
	#----------------------------- Plot the picture  -----------------------------
	require(graphics)   
	jpeg(filename=sprintf("tauA=%0.1f,tauB=%0.1f,N=%0.0f.jpeg",True_Tau_A,True_Tau_B,N),width = 672, height = 672,quality = 100) 
	qqplot(runif(200,0,1),P_chi, xlim=c(0,1),ylim=c(0,1),xlab=expression(Unif(0,1)),ylab=expression(paste("Scaled ", T[tau])))   
	title(main=substitute(paste(tau[A],"=",tauA, " ",tau[B],"=",tauB," N=",N),list(tauA=True_Tau_A,tauB=True_Tau_B,N=N)))
	a <- 0:1
	b <- a
	lines(a,b) 
	dev.off()
	return(list(biasA=biasA,varA=varA,biasB=biasB,varB=varB))
}
	
N 			<- c(100,500) 								# The number of the objects
True_Tau_A 	<- c(0.2,5)    								# True value for TauA
True_Tau_B 	<- c(0.5,10)								# True value for TauB
bA <- array(0,dim <- c(2,2))
vA <- array(0,dim <- c(2,2))
bB <- array(0,dim <- c(2,2))
vB <- array(0,dim <- c(2,2))
for (i in 1:2){
	for (j in 1:2){
		result <- Run_parameter(N[i],True_Tau_A[j],True_Tau_B[j])	
		bA[i,j] <- result$biasA
		vA[i,j] <- result$varA
		bB[i,j] <- result$biasB
		vB[i,j] <- result$varB
	}
}
sink("output.txt")
result
sink()
