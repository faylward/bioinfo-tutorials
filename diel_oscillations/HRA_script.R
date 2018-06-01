

library(snow)

################################################################################################

##### this code assumed we're looking for 24-hr periodicity. If your units are months and not hours, and you're looking for 12-month periodicity you want to change some things below (sample.times/24 -> sample.times/12, phase*24 -> phase*12, etc).
gene.glm.periodic <- function (gene.counts, Total.Counts, sample.times, minrep, maxrep) {

	Result <- list()

	gene.test <- glm(as.numeric(gene.counts) ~ offset(log(Total.Counts)) + I(cos(2*pi*sample.times/24)) + I(sin(2*pi*sample.times/24)), family=poisson)
	LR.obs <- gene.test$null.deviance - gene.test$deviance
	
	coef <- gene.test$coefficients
	amp <- sqrt( coef[2]^2 + coef[3]^2 )
	phase <- atan( -1 * coef[3] / coef[2] )

	phase <- atan( -1 * coef[3] / coef[2] )
	if (phase > 0) { peak <- 24 - phase*24/(2*pi) }
	else {  peak <- -1*phase*24/(2*pi) }
	
	if (coef[2] < 0) { 
		if (peak < 12 ) {
			peak <- peak + 12
		} else {
			peak <- peak - 12
		}
	}

	Result$fit <- gene.test
	Result$amp <- amp
	Result$peak <- peak
	Result$converged <- gene.test$converged
	Result$pval <-  pchisq(LR.obs, 2, lower.tail=FALSE)
	Result$intercept <- coef[1]
	Result$pseudoR2 <- 1- (gene.test$deviance / gene.test$null.deviance)

	perm.pass <- 0
	permutations <- 0
	perm.LR <- numeric(length=maxrep)

	while ( permutations < maxrep & (perm.pass <= 10 | permutations < minrep) ) {
		
		permutations <- permutations + 1 
		rep.times <- sample(sample.times)
		rep.test <- glm(as.numeric(gene.counts) ~ offset(log(Total.Counts)) + I(cos(2*pi*rep.times/24)) + I(sin(2*pi*rep.times/24)), family=poisson)
		rep.LR <- rep.test$null.deviance - rep.test$deviance
		perm.LR[permutations] <- rep.LR

		if ( rep.LR >= (LR.obs - 1E-5*LR.obs) )	{ perm.pass <- perm.pass+1 }

	}
	
	Result$permutations <- permutations
	Result$perm.pval <- perm.pass / permutations
	Result$perm.LR <- perm.LR[1:permutations]

	return(Result)

}

summarize.results <- function(Data.counts, Regression.Details, numgenes) {
	Result <- list()
	Result$amp <- numeric(numgenes)
	Result$peak <- numeric(numgenes)
	Result$intercept <- numeric(numgenes)
	Result$converged <- logical(numgenes)
	Result$pseudoR2 <- numeric(numgenes)
	Result$pval <- numeric(numgenes)
	Result$permutations <- numeric(numgenes)
	Result$perm.pval <- numeric(numgenes)

	for (i in 1:numgenes) {
		Result$amp[i] <- Regression.Details[[i]]$amp
		Result$peak[i] <- Regression.Details[[i]]$peak
		Result$converged[i] <- Regression.Details[[i]]$converged
		Result$pval[i] <- Regression.Details[[i]]$pval
		Result$permutations[i] <- Regression.Details[[i]]$permutations
		Result$perm.pval[i] <- Regression.Details[[i]]$perm.pval
		Result$intercept[i] <- Regression.Details[[i]]$intercept
		Result$pseudoR2[i] <- Regression.Details[[i]]$pseudoR2
	}

	Result$pval.FDR <- p.adjust(Result$pval, method="fdr")
	Result$perm.FDR <- p.adjust(Result$perm.pval, method="fdr")

	Result$table <- data.frame(Cluster=row.names(Data.Counts), Data.counts[1:numgenes,1:7], Amplitude=Result$amp, Peak.Time=Result$peak, Intercept=Result$intercept, Converged=Result$converged, PseudoR2=Result$pseudoR2, Regression.pval=Result$pval, Regression.FDR=Result$pval.FDR, Permutations=Result$permutations, Perm.pval=Result$perm.pval, Perm.FDR=Result$perm.FDR)

	print(paste(c(sum(Result$pval.FDR <= 0.1 & Result$perm.FDR <= 0.1), " out of ", numgenes, " significant."), collapse=""))

	return(Result)

}

gene.glm.getfitted <- function(Data.Counts, Regression.Details) {
	Result.table <- Data.Counts
	Total.Counts <- as.numeric(apply(Data.Counts, 2, sum))
	numgenes <- dim(Data.Counts)[1]
	for (i in 1:numgenes) {
		Result.table[i,] <- Regression.Details[[i]]$fit$fitted.values / Total.Counts
	}
	return(Result.table)
}


#######################################################################################################SAR110.25:

#### Provide a time vector that gives the times (in hours, months, or whatever as long as it's consistent) for each column of the count file to be loaded). It shouldn't matter if the spacing is uneven as long as that is reflected in the times given here. 
times=c(2, 6, 8, 12, 16, 20, 24, 28, 32, 36, 40, 48, 52, 56, 62, 66, 70, 72, 74, 78, 82, 86, 90, 94, 96, 98, 100, 104, 106, 110, 114, 118, 122, 126, 128)

# Minimum and maximum permutations to be performed for permutation tests
minrep <- 500
maxrep <- 50000

# Number of samples (columns) in count file
numsamples <- 35

# Load count file
SAR11.counts <- read.table(file="SAR406.txt", sep="\t", header=TRUE, row.names=1)
# give number of genes (or taxa, or whatever is listed as rows).
numgenes <- dim(SAR11.counts)[1]

# Get count data- this is useful if there is other metadata in the count file loaded
Data.Counts <- SAR11.counts[1:numgenes,1:35]
# Get library-specific offsets to use
Total.Counts <- as.numeric(apply(SAR11.counts[,1:35], 2, sum))

# decide how many clusters (threads) you want to run- on the server- ~6-8 is recommended. 
cl <- makeCluster(6, "SOCK")
Regression.Details <- parRapply(cl, Data.Counts, gene.glm.periodic, Total.Counts=Total.Counts, sample.times=times, minrep=minrep, maxrep=maxrep)
stopCluster(cl)
Regression.summary <- summarize.results(SAR11.counts, Regression.Details, numgenes)
write.table(Regression.summary$table, file="SAR406_HRA.results", sep="\t", row.names=FALSE, quote=FALSE)

