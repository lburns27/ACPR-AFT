#------------------------------------------------------------------------------------------------------
# R script containing algorithm for ACPR-AFT (sAFT)
#	Step 1: Pre-filter data using a marginal screening approach (call this data 'dat')
#	Steps 2 & 3: Run the Survival Prediction Algorithm (GenF & sAFT)
#	Step 4 & 5: Use results from that algorithm to estimate predicted survival (PI)
#		and calculate various measures of prediction accuracy 
#
#	Required packages: lss, flexsurv, survival, survivalROC, concreg, Hmisc
#------------------------------------------------------------------------------------------------------

library(lss)
library(flexsurv)
library(survival)
library(survivalROC)
library(concreg)
library(Hmisc)

#------------------------------------------------------------------------------------------------------
#	Steps 2 & 3: Run the Survival Prediction Algorithm  using 25 runs of training/test sets
#	Note:  The results across the 25 runs are saved as .csv files for importation later
#------------------------------------------------------------------------------------------------------

for(doit in 1:25){
	oneThird <- round(nrow(dat)/3,0)
	set.seed(doit+2000)
	trainSubj <- sample(1:nrow(dat), oneThird)
	train <- dat[-trainSubj,]
	test <- dat[trainSubj,]

	#Save training & test sets across each run
	write.csv(train, paste("Train", doit, "csv", sep = "."))
	write.csv(test, paste("Test", doit, "csv", sep = "."))

	#-------------------------------------------------------------------------------------------------
	# TRAIN SET: Loop through alpha choices (.25, .5, .75, .95) and 
	#		collect minimum PRESS & # components for each case (GF & SA)
	#-------------------------------------------------------------------------------------------------
	
	x <- as.matrix(train[,3:ncol(train)])
	y <- train$time
	K <- 15

	pressOut <- press1(x, log(y), K)				
	numComps <- which(pressOut == min(pressOut[-1]))  

	allPressK <- c()
	all_v_star_gf <- c()
	all_v_star_sa <- c()

	for(alpha in c(.25, .5, .75, .95)){

		components <- doCPR(x, log(y), alpha, a = numComps)	#Run PLS and get components

		dat_uc <- train[which(train$censor == 1),]		#Get subset of uncensored data
		x_uc <- as.matrix(dat_uc[,3:ncol(dat_uc)])
		y_uc <- dat_uc[,1]
	
		components_uc <- doCPR(x_uc, log(y_uc), alpha, a = numComps)	#Run PLS and get components

		# Fit generalized-F model and obtain model coefficients
		# Note: For any convergance issues, please see 'flexsurv' package documentation for suggestions

		modelFit_uc_gf <- flexsurvreg(Surv(time, censor) ~ components_uc, data = dat_uc, dist = "genf.orig")
		beta_uc_gf <- modelFit_uc_gf$res[,1][c(1,5:(numComps+4))]
		sigma_uc_gf <- modelFit_uc_gf$res[,1][2]
		n_uc_gf <- modelFit_uc_gf$res[,1][3]
		m_uc_gf <- modelFit_uc_gf$res[,1][4]

		#Fit semiparametric AFT model and obtain model coefficients 

		modelFit_uc_sa <- lss(Surv(log(time), censor) ~ components_uc, data = dat_uc, gehanonly = T, mcsize = 200, maxiter = 50)
		beta_uc_sa <- modelFit_uc_sa$betag

		v <- log(train$time)

		#GF adjustement ------------------------------------------------------------------------------------
		comp2 <- cbind(rep(1,nrow(train)), components)	
		v_star_gf <- rep(0, nrow(train))

		for(i in 1:nrow(train)){
			if(train$censor[i] == 1){v_star_gf[i] <- v[i]}					
			else{	v_star_gf[i] <- beta_uc_gf%*%comp2[i,] + sigma_uc_gf*expectGF(i)	
			}
		}	

		#SA adjustment ------------------------------------------------------------------------------------
		e <- v - components%*%beta_uc_sa				
		eps <- .Machine$double.eps^(2/3)			
		zdummy <- matrix(rep(1, length(v)), ncol = 1)	
		eres <- lss.eres(e, zdummy, train$censor, eps)	
		v_star_sa <- v*train$censor + (1 - train$censor)*(eres + components%*%beta_uc_sa)


		#PRESS & K based on adjusted times ------------------------------------------------------------------------------------

		v_star_gf <- as.vector(v_star_gf)
		v_star_sa <- as.vector(v_star_sa)

		pressOut_gf <- press1(x, v_star_gf, K)				
		numComps_gf<- which(pressOut_gf == min(pressOut_gf[-1])) 
		minPress_gf<- min(pressOut_gf[-1])

		pressOut_sa <- press1(x, v_star_sa, K)				
		numComps_sa <- which(pressOut_sa == min(pressOut_sa[-1])) 
		minPress_sa<- min(pressOut_sa[-1])


		allPressK <- rbind(allPressK, c(numComps, alpha, minPress_gf, numComps_gf, minPress_sa, numComps_sa))
		all_v_star_gf <- cbind(all_v_star_gf, v_star_gf)
		all_v_star_sa <- cbind(all_v_star_sa, v_star_sa)

		#Save results for each run
		write.csv(allPressK, paste("AllPRESS", doit, "csv", sep = "."))    
		write.csv(all_v_star_gf, paste("AllVstarGF", doit, "csv", sep = "."))    
		write.csv(all_v_star_sa, paste("AllVstarSA", doit, "csv", sep = "."))    
		}
	}


#------------------------------------------------------------------------------------------------------
#	Steps 4 & 5: Use results from that algorithm to estimate predicted survival (PI)
#		and calculate various measures of prediction accuracy
#------------------------------------------------------------------------------------------------------

#Specify ttROC = 2, 3 & 5 year time points
# Example: ttROC <- c(2,3,5) if measured in days; ttROC <- c(24, 36,60) if measured in months

#Initialize results 
GFresults <- c()
SAresults <- c()
UnadjGFresults <- c()
UnadjSAresults <- c()

for(doit in 1:25){
	# Import results from 25 runs ---------------------------------------------------------------------
	alphaK <- read.csv(paste("AllPRESS", doit, "csv", sep = "."), row.names = 1)
	v_star_gf <- read.csv(paste("AllVstarGF", doit, "csv", sep = "."), row.names = 1)
	v_star_sa <- read.csv(paste("AllVstarSA", doit, "csv", sep = "."), row.names = 1)
	train <- read.csv(paste("Train", doit, "csv", sep = "."), row.names = 1)
	test <- read.csv(paste("Test", doit, "csv", sep = "."), row.names = 1)

	#Select optimal (alpha, k) combo---------------------------------------------------------------------

	numComps <- alphaK[1,1]
	minGF <- which(alphaK[,3] == min(alphaK[,3]))
	alpha_gf <- alphaK[,2][minGF]
	numComps_gf <- alphaK[,4][minGF]
	minSA <- which(alphaK[,5] == min(alphaK[,5]))
	alpha_sa <- alphaK[,2][minSA]
	numComps_sa <- alphaK[,6][minSA]
	v_star_gf <- as.vector(v_star_gf[,minGF])
	v_star_sa <- as.vector(v_star_sa[,minSA])


	# Get PLS coeff & analyze train set results-------------------------------------------------------------

	x <- as.matrix(train[,3:ncol(train)])
	coeff_gf <- doCPR(x, v_star_gf, alpha_gf, a = numComps_gf, LOOCV = TRUE)
	coeff_sa <- doCPR(x, v_star_sa, alpha_sa, a = numComps_sa, LOOCV = TRUE)

	#Calculate PI for training set & test set

	pred_train_gf <- x%*%coeff_gf
	pred_train_sa <- x%*%coeff_sa

	x_test <- as.matrix(test[,3:ncol(test)])
	pred_test_gf <- x_test%*%coeff_gf
	pred_test_sa <- x_test%*%coeff_sa

	#Get R^2 & MSE---------------------------------------------
	R2_gf <- doCPR(x, v_star_gf, alpha_gf, a = numComps_gf, R2out = TRUE)
	R2_sa <- doCPR(x, v_star_sa, alpha_sa, a = numComps_sa, R2out = TRUE)
	MSE_train_gf <- sum(train$censor*(pred_train_gf - log(train$time))^2)/sum(train$censor)
	MSE_train_sa <- sum(train$censor*(pred_train_sa - log(train$time))^2)/sum(train$censor)
	MSE_test_gf <- sum(test$censor*(pred_test_gf - log(test$time))^2)/sum(test$censor)
	MSE_test_sa <- sum(test$censor*(pred_test_sa - log(test$time))^2)/sum(test$censor)
	

	#ROC curves---------------------------------------------------------------------------------

	auc_gf2 <- survivalROC(test$time, test$censor, marker = pred_test_gf, predict.time = ttROC[1], method = "KM")$AUC
	auc_gf3 <- survivalROC(test$time, test$censor, marker = pred_test_gf, predict.time = ttROC[2], method = "KM")$AUC
	auc_gf5 <- survivalROC(test$time, test$censor, marker = pred_test_gf, predict.time = ttROC[3], method = "KM")$AUC

	auc_sa2 <- survivalROC(test$time, test$censor, marker = pred_test_sa, predict.time = ttROC[1], method = "KM")$AUC
	auc_sa3 <- survivalROC(test$time, test$censor, marker = pred_test_sa, predict.time = ttROC[2], method = "KM")$AUC
	auc_sa5 <- survivalROC(test$time, test$censor, marker = pred_test_sa, predict.time = ttROC[3], method = "KM")$AUC

	auc_gf2TR <- survivalROC(train$time, train$censor, marker = pred_train_gf, predict.time = ttROC[1], method = "KM")$AUC
	auc_gf3TR <- survivalROC(train$time, train$censor, marker = pred_train_gf, predict.time = ttROC[2], method = "KM")$AUC
	auc_gf5TR <- survivalROC(train$time, train$censor, marker = pred_train_gf, predict.time = ttROC[3], method = "KM")$AUC

	auc_sa2TR <- survivalROC(train$time, train$censor, marker = pred_train_sa, predict.time = ttROC[1], method = "KM")$AUC
	auc_sa3TR <- survivalROC(train$time, train$censor, marker = pred_train_sa, predict.time = ttROC[2], method = "KM")$AUC
	auc_sa5TR <- survivalROC(train$time, train$censor, marker = pred_train_sa, predict.time = ttROC[3], method = "KM")$AUC

	#Concreg & Harrel's C

	concMod_TR_gf <- concreg(data = dat, Surv(dat$time, censor) ~ pred_train_gf)
	concTRgf <- .5 + abs(cindex(concMod_TR_gf)-.5)
	concMod_TR_sa <- concreg(data = dat, Surv(dat$time, censor) ~ pred_train_sa)
	concTRsa <- .5 + abs(cindex(concMod_TR_sa)-.5)
	concMod_TE_gf <- concreg(data = dat, Surv(dat$time, censor) ~ pred_test_gf)
	concTEgf <- .5 + abs(cindex(concMod_TE_gf)-.5)
	concMod_TE_sa <- concreg(data = dat, Surv(dat$time, censor) ~ pred_test_sa)
	concTEsa <- .5 + abs(cindex(concMod_TE_sa)-.5)

	cTRgf <- rcorrcens(Surv(dat$time, dat$censor) ~ pred_train_gf)[1]
	cTRsa <- rcorrcens(Surv(dat$time, dat$censor) ~ pred_train_sa)[1]
	cTEgf <- rcorrcens(Surv(dat$time, dat$censor) ~ pred_test_gf)[1]
	cTEsa <- rcorrcens(Surv(dat$time, dat$censor) ~ pred_test_sa)[1]

	# Unadjusted analysis ----------------------------------------------------------------------------

	y <- train$time
	R2_unadj_gf <- doCPR(x, log(y), alpha_gf, a = numComps, R2out = TRUE)
	R2_unadj_sa <- doCPR(x, log(y), alpha_sa, a = numComps, R2out = TRUE)

	coeff_unadj_gf <- doCPR(x, log(y), alpha_gf, a = numComps, LOOCV = TRUE)
	coeff_unadj_sa <- doCPR(x, log(y), alpha_sa, a = numComps, LOOCV = TRUE)

	pred_train_unadj_gf <- x%*%coeff_unadj_gf
	pred_train_unadj_sa <- x%*%coeff_unadj_sa
	MSE_train_unadj_gf <- sum(train$censor*(pred_train_unadj_gf - log(train$time))^2)/sum(train$censor)
	MSE_train_unadj_sa <- sum(train$censor*(pred_train_unadj_sa - log(train$time))^2)/sum(train$censor)

	pred_test_unadj_gf <- x_test%*%coeff_unadj_gf
	pred_test_unadj_sa <- x_test%*%coeff_unadj_sa
	MSE_test_unadj_gf <- sum(test$censor*(pred_test_unadj_gf - log(test$time))^2)/sum(test$censor)
	MSE_test_unadj_sa <- sum(test$censor*(pred_test_unadj_sa - log(test$time))^2)/sum(test$censor)

	auc_unadj_gf2 <- survivalROC(test$time, test$censor, marker = pred_test_unadj_gf, predict.time = ttROC[1], method = "KM")$AUC
	auc_unadj_gf3 <- survivalROC(test$time, test$censor, marker = pred_test_unadj_gf, predict.time = ttROC[2], method = "KM")$AUC
	auc_unadj_gf5 <- survivalROC(test$time, test$censor, marker = pred_test_unadj_gf, predict.time = ttROC[3], method = "KM")$AUC

	auc_unadj_sa2 <- survivalROC(test$time, test$censor, marker = pred_test_unadj_sa, predict.time = ttROC[1], method = "KM")$AUC
	auc_unadj_sa3 <- survivalROC(test$time, test$censor, marker = pred_test_unadj_sa, predict.time = ttROC[2], method = "KM")$AUC
	auc_unadj_sa5 <- survivalROC(test$time, test$censor, marker = pred_test_unadj_sa, predict.time = ttROC[3], method = "KM")$AUC

	auc_unadj_gf2TR <- survivalROC(train$time, train$censor, marker = pred_train_unadj_gf, predict.time = ttROC[1], method = "KM")$AUC
	auc_unadj_gf3TR <- survivalROC(train$time, train$censor, marker = pred_train_unadj_gf, predict.time = ttROC[2], method = "KM")$AUC
	auc_unadj_gf5TR <- survivalROC(train$time, train$censor, marker = pred_train_unadj_gf, predict.time = ttROC[3], method = "KM")$AUC

	auc_unadj_sa2TR <- survivalROC(train$time, train$censor, marker = pred_train_unadj_sa, predict.time = ttROC[1], method = "KM")$AUC
	auc_unadj_sa3TR <- survivalROC(train$time, train$censor, marker = pred_train_unadj_sa, predict.time = ttROC[2], method = "KM")$AUC
	auc_unadj_sa5TR <- survivalROC(train$time, train$censor, marker = pred_train_unadj_sa, predict.time = ttROC[3], method = "KM")$AUC

	#Concreg & Harrel's C
	concMod_TR_gf_unadj <- concreg(data = dat, Surv(dat$time, censor) ~ pred_train_unadj_gf)
	concTRgf_unadj <- .5 + abs(cindex(concMod_TR_gf_unadj)-.5)
	concMod_TR_sa_unadj <- concreg(data = dat, Surv(dat$time, censor) ~ pred_train_unadj_sa)
	concTRsa_unadj <- .5 + abs(cindex(concMod_TR_sa_unadj)-.5)
	concMod_TE_gf_unadj <- concreg(data = dat, Surv(dat$time, censor) ~ pred_test_unadj_gf)
	concTEgf_unadj <- .5 + abs(cindex(concMod_TE_gf_unadj)-.5)
	concMod_TE_sa_unadj <- concreg(data = dat, Surv(dat$time, censor) ~ pred_test_unadj_sa)
	concTEsa_unadj <- .5 + abs(cindex(concMod_TE_sa_unadj)-.5)

	cTRgf_unadj <- rcorrcens(Surv(dat$time, dat$censor) ~ pred_train_unadj_gf)[1]
	cTRsa_unadj <- rcorrcens(Surv(dat$time, dat$censor) ~ pred_train_unadj_sa)[1]
	cTEgf_unadj <- rcorrcens(Surv(dat$time, dat$censor) ~ pred_test_unadj_gf)[1]
	cTEsa_unadj <- rcorrcens(Surv(dat$time, dat$censor) ~ pred_test_unadj_sa)[1]

	# Summarize all results across runs  -----------------------------------------------------------------------		
	
	GFresults <- rbind(GFresults, c(R2_gf[numComps_gf], MSE_train_gf, auc_gf2TR, auc_gf3TR, auc_gf5TR, concTRgf, cTRgf, MSE_test_gf, auc_gf2, auc_gf3, auc_gf5, concTEgf, cTEgf))
	SAresults <- rbind(SAresults, c(R2_sa[numComps_sa], MSE_train_sa, auc_sa2TR, auc_sa3TR, auc_sa5TR, concTRsa, cTRsa, MSE_test_sa, auc_sa2, auc_sa3, auc_sa5, concTEsa, cTEsa))
	UnadjGFresults <- rbind(UnadjGFresults, c(R2_unadj_gf[numComps], MSE_train_unadj_gf, auc_unadj_gf2TR, auc_unadj_gf3TR, auc_unadj_gf5TR, concTRgf_unadj, cTRgf_unadj,
						MSE_test_unadj_gf, auc_unadj_gf2, auc_unadj_gf3, auc_unadj_gf5, concTEgf_unadj, cTEgf_unadj))
	UnadjSAresults <- rbind(UnadjSAresults, c(R2_unadj_sa[numComps], MSE_train_unadj_sa, auc_unadj_sa2TR, auc_unadj_sa3TR, auc_unadj_sa5TR, concTRsa_unadj, cTRsa_unadj,
						MSE_test_unadj_sa, auc_unadj_sa2, auc_unadj_sa3, auc_unadj_sa5, concTEsa_unadj, cTEsa_unadj))

}

#Calculate median results
apply(GFresults, 2, function(x) median(x, na.rm = T))
apply(SAresults, 2, function(x) median(x, na.rm = T))
apply(UnadjGFresults, 2, function(x) median(x, na.rm = T))
apply(UnadjSAresults, 2, function(x) median(x, na.rm = T))

