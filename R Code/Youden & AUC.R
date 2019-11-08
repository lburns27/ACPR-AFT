#----------------------------------------------------------------------------
# getVIPandCoeffResults: Function to calculate Sens, Spec, Youden & AUC
#	for VIP & PLS coeff
#	
#	out: Data frame with n rows (number of genes) and 2 columns ($cprCoeff and $vip)
#	effectGenes: number of significant genes 
#
#	Returns: Sensitivity, Specificity, Youden, AUC & AUC standard deviation
#		for CPR coefficients & VIP
#
#	Required package: MESS
#-----------------------------------------------------------------------------

library(MESS)

getVIPandCoeffResults <- function(out, effectGenes){
	n <- effectGenes
	m <- nrow(out)

	goodgenes <- rownames(out)[1:n]
	badgenes <- rownames(out)[(n+1):nrow(out)]

	sortBeta <- out[order(abs(out$cprCoeff), decreasing = T),]
	sortVIP <- out[order(abs(out$vip), decreasing = T),]

	#Sens & Spec
	sb <- rownames(sortBeta)[1:n]
	sv <- rownames(sortVIP)[1:n]

	list1 <- list(sb, sv)
	tpr <- c()
	fpr <- c()
	for(j in 1:length(list1)){
		int <- list1[[j]]
		goodpick <- length(int[! int %in% badgenes])
		badpick <- length(int[int %in% badgenes])
		tpr <- c(tpr, goodpick/n)
		fpr <- c(fpr, 1-badpick/(m-n))
		}

	both <- as.data.frame(cbind(tpr, fpr, tpr+fpr-1))
	colnames(both) <- c("Sensitivity (TPR)", "Specificity (TNR)", "Youden")
	rownames(both) <- c("PLS_Coeff", "VIP")


	choosenum <- seq(0, m, 10)
	tpr <- matrix(0, length(choosenum), 2) 
	fpr <- matrix(0, length(choosenum), 2) 
	colnames(tpr) <- colnames(fpr) <- c("beta", "vip")

	count <- 0
	for(i in choosenum){
		count <- count+1
		b <- rownames(sortBeta)[1:i]
		v <- rownames(sortVIP)[1:i]

		list1 <- list(b,v)

		for(j in 1:length(list1)){
			int <- list1[[j]]
			goodpick <- length(int[! int %in% badgenes])
			badpick <- length(int[int %in% badgenes])
			tpr[count,j] <- goodpick/n
			fpr[count,j] <- badpick/(m-n)
		}
	}


	area <- c()

	for(i in 1:2){
		area <- c(area, auc(fpr[,i], tpr[,i], 0, 1) )
		}

	changeArea <- function(area, lower, upper){
		area[1] <- auc(fpr[,1], tpr[,1], lower, upper)
		area[2] <- auc(fpr[,2], tpr[,2], lower, upper)
		return(area)
	}
	
	if(sum(is.na(area)) > 0){area <- changeArea(area, .001, 1)}
	if(sum(is.na(area)) > 0){area <- changeArea(area, .01, .9999)}
	if(sum(is.na(area)) > 0){area <- changeArea(area, .01, .999)}
	if(sum(is.na(area)) > 0){area <- changeArea(area, .01, .99)}
	if(sum(is.na(area)) > 0){area <- changeArea(area, .01, .9)}


	stdevAUC <- function(x){
		q1 <- x/(2-x)
		q2 <- 2*x^2/(1+x)
		se <- sqrt((x*(1-x)+(n-1)*(q1-x^2)+(m-n-1)*(q2-x^2))/(n*(m-n)))
		return(se)
		}

	sdAUC <- stdevAUC(area)
	
	both$AUC <- area
	both$AUC_stdev <- sdAUC
	
	spitOut <- matrix(c(both[1,], both[2,]), 1, 10)

	return(spitOut)
}

