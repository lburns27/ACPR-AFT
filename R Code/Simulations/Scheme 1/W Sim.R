#----------------------------------------------------------------------------------------
# Data generating function: W (Weibull, same shape, different scale)
#	N: Number of subjects
#	small: # of genes with small effect
#	large: # of genes with large effect
#	tau: censoring parameter (see examples below)
#	NoSig: number of genes with no effect
#	lam1, lam2: scale parameters
#----------------------------------------------------------------------------------------


makeDataW <- function(N, small, large, tau, NoSig, lam1, lam2){
	Bo_small <- .599
	Bo_large <- 2

	v <- runif(N, min=0, max=1)
	y <- (-log(v))^(1/2)/0.5

	if(tau != "inf"){
		z <- runif(N, 0.00001, tau)
		t <- pmin(y, z)
		d <- (y <= z) + 0
		}
	if(tau == "inf"){ 
		t <- y
		d <- rep(1,N) 
		}

	survDat_unsorted <- cbind(t, d)
	survDat <- as.data.frame(survDat_unsorted[order(t),])
	dat <- survDat

	Bo <- Bo_small
	Bt <- function(time) {rep(Bo*log(lam1/lam2), length(time))} 
	for(i in 1:small){
		a <- rnorm(N, 0, 1)
		x <- c()
		R <- 1:N

		for(i in 1:N){
			if(length(R) == 1){ x[i] <- a[R] }
			else if (survDat$d[i] == 0){
				j <- sample(R, 1)
				x[i] <- a[j]
				R <- R[R != j]
				}
			else{
				probs <- exp(Bt(survDat$t[i])*a[R])/sum(exp(Bt(survDat$t[i])*a[R]))
				j <- sample(R, 1, prob = probs)
				x[i] <- a[j]
				R <- R[R != j]
				}
			}
		dat <- cbind(dat, x)
		}

	Bo <- Bo_large
	Bt <- function(time) {rep(Bo*log(lam1/lam2), length(time))} 
	for(i in 1:large){
		a <- rnorm(N, 0, 1)
		x <- c()
		R <- 1:N

		for(i in 1:N){
			if(length(R) == 1){ x[i] <- a[R] }
			else if (survDat$d[i] == 0){
				j <- sample(R, 1)
				x[i] <- a[j]
				R <- R[R != j]
				}
			else{
				probs <- exp(Bt(survDat$t[i])*a[R])/sum(exp(Bt(survDat$t[i])*a[R]))
				j <- sample(R, 1, prob = probs)
				x[i] <- a[j]
				R <- R[R != j]
				}
			}
		dat <- cbind(dat, x)
		}

	Bo <- 0
	Bt <- function(time) {rep(Bo*log(lam1/lam2), length(time))} 
	for(i in 1:NoSig){
		a <- rnorm(N, 0, 1)
		x <- c()
		R <- 1:N

		for(i in 1:N){
			if(length(R) == 1){ x[i] <- a[R] }
			else if (survDat$d[i] == 0){
				j <- sample(R, 1)
				x[i] <- a[j]
				R <- R[R != j]
				}
			else{
				probs <- exp(Bt(survDat$t[i])*a[R])/sum(exp(Bt(survDat$t[i])*a[R]))
				j <- sample(R, 1, prob = probs)
				x[i] <- a[j]
				R <- R[R != j]
				}
			}
		dat <- cbind(dat, x)
		}


	colnames(dat) <- c("time", "censor", paste("S", 1:small, sep = ""), 
				paste("L", 1:large, sep = ""), paste("NoSig", 1:NoSig, sep = ""))
	return(dat)	
}


# Example----------------------------------------------------------------------------------------------------------

numSmall <- 200
numLarge <- 200
notsig <- 4600
N <- 200

# 0% censoring
dat0 <- makeDataW(N, numSmall, numLarge, "inf", notsig, 1,.5)

# 33% censoring
dat33 <- makeDataW(N, numSmall, numLarge, 5.3, notsig, 1,.5)

# 67% censoring
dat67 <- makeDataW(N, numSmall, numLarge, 2.4, notsig, 1,.5)

# 80% censoring
dat80 <- makeDataW(N, numSmall, numLarge, 1.7, notsig, 1,.5)

