#----------------------------------------------------------------------------------
# R script containing algorithm for CPR-AFT
#	Step 1: Import data & define x & y (x = gene matrix, y = survival time)
#	Step 2: Get functions for CPR & LOOCV based on PRESS
#	*Note: Found in the "CPR & LOOCV based on PRESS" R file in our Github repository
#	Step 3: Run CPR-AFT
#----------------------------------------------------------------------------------

	x <- as.matrix(dat[,3:ncol(dat)])
	y <- dat$time

# Run LOOCV based on PRESS to find optimal number of components (from 2 to 15)
# Note: K = 15 could be changed if desired to an alternate maximum value

	K <- 15

	pressOut <- press1(x, log(y), K)				
	numComps <- which(pressOut == min(pressOut[-1]))  

# Get components and CPR coefficients based on optimal numComps chosen above

	#Any alpha value can be chosen
	alpha <- 0.5 

	components <- doCPR(x, log(y), alpha, a = numComps)
	beta <- doCPR(x, log(y), alpha, a = numComps, beta = TRUE)

# Options for further work:
#	1. The components can be utilized for additional analyses.
#	2. The CPR coefficients (beta) can be used to predict survival time (see below).
#	3. Get VIP for each gene to rank genes (see below) 		

# Option 2
	PI <- x%*%beta 

# Option 3
	vip <- doCPR(x, log(y), alpha, a = numComps, VIP = TRUE)






	
