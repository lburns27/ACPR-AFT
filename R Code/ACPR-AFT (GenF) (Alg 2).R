#----------------------------------------------------------------------------------
# R script containing algorithm for ACPR-AFT (GenF)
#	Step 1: Import data (named 'dat') & define x & y (x = gene matrix, y = survival time)
#	Step 2: Get functions for CPR & LOOCV based on PRESS
#	*Note: Found in the "CPR & LOOCV based on PRESS" R file in our Github repository
#	Step 3: Run ACPR-AFT (GenF)
#
#	Required packages: flexsurv
#----------------------------------------------------------------------------------

library(flexsurv)

# Run expectGF function:  Used to compute adjusted times later in step 3 of algorithm 
	
	expectGF <- function(whatObs){
		mu_uc_gf <- beta_uc_gf%*%comp2[whatObs,]
		val <- (v[whatObs] - mu_uc_gf) / sigma_uc_gf
		int_lim <- (n_uc_gf/m_uc_gf)*exp(val)/(1 + (n_uc_gf/m_uc_gf)*exp(val))
			
		betafunc <- function(x){ factorial(n_uc_gf + m_uc_gf-1)/(factorial(n_uc_gf-1)*factorial(m_uc_gf-1)) * x^(n_uc_gf-1)*(1-x)^(m_uc_gf-1) }
		betafunc2 <- function(x){ log((m_uc_gf/n_uc_gf)*x/(1-x))*
				factorial(n_uc_gf + m_uc_gf -1)/(factorial(n_uc_gf-1)*factorial(m_uc_gf-1))*x^(n_uc_gf-1)*(1-x)^(m_uc_gf-1)}

		numerator <- integrate(betafunc2, int_lim, 1)$val
		denominator <- 1 - integrate(betafunc, 0, int_lim)$val
		expect <- numerator/denominator
		return(expect)
	}

# Define x & y

	x <- as.matrix(dat[,3:ncol(dat)])
	y <- dat$time

# STEP 1:  Obtain initial optimal number of components using LOOCV based on PRESS
# Note: K = 15 could be changed if desired to an alternate maximum value

	K <- 15
	pressOut <- press1(x, log(y), K)				
	numComps <- which(pressOut == min(pressOut[-1]))  

# STEPS 2-5: Run ACPR-AFT for various alpha values
# 	Option 1: Run for a single alpha value (choose 'aa' as single value)
#	Option 2: Run for a selection of alpha values (.25, .5, .75, .95) & choose 
#		optimal alpha based on lowest PRESS; Choose aa <- c(.25, .5, .75, .95)

	aa <- c(.25, .5, .75, .95)

	for(alpha in aa){

		components <- doCPR(x, log(y), alpha, a = numComps)	#Run CPR and get components

		dat_uc <- dat[which(dat$censor == 1),]			#Get subset of uncensored data
		x_uc <- as.matrix(dat_uc[,3:ncol(dat_uc)])
		y_uc <- dat_uc[,1]
	
		components_uc <- doCPR(x_uc, log(y_uc), alpha, a = numComps)	#Run CPR and get components

		# Fit generalized-F model and obtain model coefficients
		# Note: For any convergance issues, please see 'flexsurv' package documentation for suggestions

		modelFit_uc_gf <- flexsurvreg(Surv(time, censor) ~ components_uc, data = dat_uc, dist = "genf.orig")
		beta_uc_gf <- modelFit_uc_gf$res[,1][c(1,5:(numComps+4))]
		sigma_uc_gf <- modelFit_uc_gf$res[,1][2]
		n_uc_gf <- modelFit_uc_gf$res[,1][3]
		m_uc_gf <- modelFit_uc_gf$res[,1][4]

		# Adjusted survival times based on GenF model fit above

		v <- log(dat$time)

			comp2 <- cbind(rep(1,nrow(dat)), components)	
			v_star_gf <- rep(0, nrow(dat))

			for(i in 1:nrow(dat)){
				if(dat$censor[i] == 1){v_star_gf[i] <- v[i]}					
				else{	v_star_gf[i] <- beta_uc_gf%*%comp2[i,] + sigma_uc_gf*expectGF(i)	
				}
			}	

		# Run press to get optimal number of components

		v_star_gf <- as.vector(v_star_gf)
		pressOut_gf <- press1(x, v_star_gf, K)				
		numComps_gf<- which(pressOut_gf == min(pressOut_gf[-1])) 
		minPress_gf<- min(pressOut_gf[-1])

		#Store results across each alpha value 

		allPressK <- rbind(allPressK, c(numComps, alpha, minPress_gf, numComps_gf))
		all_v_star_gf <- cbind(all_v_star_gf, v_star_gf)

	}
	}

# STEP 6: Identify optimal (alpha, K) combo and compute CPR coefficients 

	minGF <- which(allPressK[,3] == min(allPressK[,3]))
	alpha_gf <- allPressK[,2][minGF]
	numComps_gf <- allPressK[,4][minGF]
	v_star_gf <- as.vector(all_v_star_gf[,minGF])

	beta_gf <- doCPR(x, v_star_gf, alpha_gf, a = numComps_gf, beta = TRUE)

	#Use the CPR coefficients to computed predicted survival times
	PI <- x%*%beta_gf

# Options for further work: 
#	1. Use the CPR components (see below) for further analyses
#	2. Get VIP for each gene to rank genes (see below) 		

# Option 1
	components_gf <- doCPR(x, v_star_gf, alpha_gf, a = numComps_gf)

# Option 2
	vip <- doCPR(x, v_star_gf, alpha_gf, a = numComps_gf, VIP = TRUE)
	

