C#----------------------------------------------------------------------------------
# R script containing functions for CPR & LOOCV based on PRESS
#	CPR function: doCPR
#	LOOCV function: press1
#----------------------------------------------------------------------------------

#----------------------------------------------------------------------------------
# doCPR: function for running CPR given alpha (a)
#	x: x matrix from data (predictors, gene matrix)
#	y: outcome (survival time)
#	alpha: specified alpha value (0.25, 0.50, 0.75, 0.95) 
#		*Note: Other alpha values can be specified in the interval (0,1)
#	a: number of components
#	LOOCV: Option TRUE returns CPR coefficients
#	R2out: Option TRUE returns R-squared
#	VIP: Option TRUE calculates and returns VIP values
#	Note: When beta = R2out = VIP = FALSE, the CPR components are returned 
#-----------------------------------------------------------------------------------

doCPR <- function(x, y, alpha, a, beta = FALSE, R2out = FALSE, VIP = FALSE){

	svd <- svd(x)
	l <- svd$d^2
	u <- svd$u
	v <- svd$v

	r <- sum(l > l[1]/1e14)

	l <- l[1:r]
	u <- u[,1:r]
	v <- v[,1:r]

	a <- min(a,r)

	gam <- alpha/(1-alpha)

	lgam <- l^gam

	rho <- t(y%*%u)

	T <- matrix(0, r, a)

	for(i in 1:a){
		tt <- lgam*rho
		tt <- tt-T[,1:max(1,i-1)]%*%(t(T[,1:max(1,i-1)])%*%tt)
		tt <- tt/sqrt(t(tt)%*%tt)[1]
		T[,i] <- tt
		rho <- rho - tt%*%(t(tt)%*%rho)
	}

	c_loadings <- t(rho)%*%T	#Doesn't do what its supposed to do....

	R <- v%*%(diag(1/sqrt(l))%*%T)	#X block weights
	P <- v%*%(diag(sqrt(l))%*%T)		#X block loadings

	#Doesn't work bc c_loadings not integers --------------------------
	B_CPR = R%*%upper.tri(t(c_loadings)*matrix(1,a))  

	#final_T <- u%*%T
	final_T <- x%*%R

	#Predict y--NOT in MATLAB code
	c2 <- t(y)%*%u%*%T
	y_predict <- x%*%R%*%t(c2)
	bpls <- R%*%t(c2)

	R2X <- 100*t(apply(T^2, 1, cumsum))
	R2y <- 100*cumsum(c2^2)/(t(y)%*%y) #Use c_loadings or c2?

	if(VIP){
		vip <- c()
		for(i in 1:nrow(R)){
   			SS <- c2^2 * colSums(final_T^2)
  			Rnorm <- colSums(R^2)
    			vip <- c(vip, sqrt(nrow(R)*sum(SS*R[i,]^2/Rnorm)/sum(SS)))
			}
		return(vip)
		}


	if(beta){ return(bpls) }
	if(R2out){return(R2y)}
	else{ return(final_T) }
	}

#----------------------------------------------------------------------------------
# press1: function for running LOOCV based on press
#	x: x matrix from data (predictors, gene matrix)
#	y: outcome (survival time)
#	a: number of maximum components to examine
#-----------------------------------------------------------------------------------

press1 <- function(x, y, numComps){
	ssdev <- rep(0, numComps)
	for(k in 2:numComps){
		for(i in 1:nrow(x)){
			x2 <- x[-i,]
			y2 <- y[-i]
			bpls <- doCPR(x2, y2, 0.5, k, beta = TRUE)
			pred <- x[i,]%*%bpls
			ssdev[k] <- ssdev[k] + (pred - y[i])^2
		}
	}
ssdev
}


