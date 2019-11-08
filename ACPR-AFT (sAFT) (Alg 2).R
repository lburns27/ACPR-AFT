#----------------------------------------------------------------------------------
# R script containing algorithm for ACPR-AFT (sAFT)
#	Step 1: Import data (named 'dat') & define x & y (x = gene matrix, y = survival time)
#	Step 2: Get functions for CPR & LOOCV based on PRESS
#	*Note: Found in the "CPR & LOOCV based on PRESS" R file in our Github repository
#	Step 3: Run ACPR-AFT (sAFT)
#
#	Required packages: lss
#----------------------------------------------------------------------------------

library(lss)

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

		#Fit semiparametric AFT model and obtain model coefficients 

		modelFit_uc_sa <- lss(Surv(log(time), censor) ~ components_uc, data = dat_uc, gehanonly = T, mcsize = 200, maxiter = 50)
		beta_uc_sa <- modelFit_uc_sa$betag

		# Adjusted survival times based on semiparametric AFT model fit above

		v <- log(dat$time)

			e <- v - components%*%beta_uc_sa				
			eps <- .Machine$double.eps^(2/3)			
			zdummy <- matrix(rep(1, length(v)), ncol = 1)	
			eres <- lss.eres(e, zdummy, dat$censor, eps)	
			v_star_sa <- v*dat$censor + (1 - dat$censor)*(eres + components%*%beta_uc_sa)

		# Run press to get optimal number of components

		v_star_sa <- as.vector(v_star_sa)
		pressOut_sa <- press1(x, v_star_sa, K)				
		numComps_sa <- which(pressOut_sa == min(pressOut_sa[-1])) 
		minPress_sa<- min(pressOut_sa[-1])

		#Store results across each alpha value 

		allPressK <- rbind(allPressK, c(numComps, alpha, minPress_sa, numComps_sa))
		all_v_star_sa <- cbind(all_v_star_sa, v_star_sa)
	}
	}

# STEP 6: Identify optimal (alpha, K) combo and compute CPR coefficients 

	minSA <- which(allPressK[,3] == min(allPressK[,3]))
	alpha_sa <- allPressK[,2][minSA]
	numComps_sa <- allPressK[,4][minSA]
	v_star_sa <- as.vector(all_v_star_sa[,minSA])

	beta_sa <- doCPR(x, v_star_sa, alpha_sa, a = numComps_sa, beta = TRUE)

	#Use the CPR coefficients to computed predicted survival times
	PI <- x%*%beta_sa

# Options for further work: 
#	1. Use the CPR components (see below) for further analyses
#	2. Get VIP for each gene to rank genes (see below) 		

# Option 1
	components_sa <- doCPR(x, v_star_sa, alpha_sa, a = numComps_sa)

# Option 2
	vip <- doCPR(x, v_star_sa, alpha_sa, a = numComps_sa, VIP = TRUE)

