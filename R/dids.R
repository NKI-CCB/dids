#data: matrix of data; genes per row, samples per column
#cl: class vector; 1 = controls/responder/sensitive; 0 = cases/non-responder/resistant
#modF: scoring function parameter. See cutf for details.
#modType: type of scoring function to use. See cutf for details.
#alt: specify one-sided or two-sided test (i.e. 'greater', 'less' or 'two.sided'). A two sided test will perform 'less' and 'greater' separately and return the highest score between the two resulting scores for each feature/gene.

DIDSscore <- function(data, cl, modF=2, modType='tanh', alt='two.sided'){
	if(alt == 'two.sided'){ return(DIDSscore_combined(data, cl, modF=modF, modType=modType)) }
	else{ return(DIDSscore_directed(data, cl, alt=alt, modF=modF, modType=modType)) }
}

#convenience function to calculate two sided DIDS scores. Use 'DIDSscore' as a wrapper.
DIDSscore_combined <- function(data, cl, modF=2, modType='tanh'){
	DIDSup <- DIDSscore(data, cl, alt='greater', modF=modF, modType=modType)
	DIDSdown <- DIDSscore(data, cl, alt='less', modF=modF, modType=modType)

	DIDScombined <- apply(cbind(DIDSup, DIDSdown),1,function(x){
		maxIndex <- which.max(abs(x[c(1,3)]))
		#silly transform to get left or right pair of statistic + p-value
		returnIndex <- (maxIndex - 1) * 2 + 1
		return(x[c(returnIndex,(returnIndex+1))])
		})

	result <- as.data.frame(t(DIDScombined))

	return(result)
}

#function to calculate actual one-sided DIDS score and corresponding p-value. Use 'DIDSscore' as a wrapper.
DIDSscore_directed <- function(data, cl, alt='greater', modF=2, modType='tanh'){
	pastThreshold <- pastThreshold(data, cl, alt)

	numberAbove <- rowSums(!is.na(pastThreshold))
	pvals <- test.tabulation(sum(cl==1), sum(cl==0))[numberAbove+1]
	pastThresholdCor <- t(apply(pastThreshold, 1, function(x)cutf(x, modF, modType)))

	score <- rowSums(pastThresholdCor, na.rm=T)

	if(alt == 'less'){ score <- -score }

	return(data.frame(score=score, pval=pvals))
}

#function to calculate the amount by which a feature is above (or below) the maximum (or minimum) of the control group
#used by DIDSscore_directed and explainedCases
pastThreshold <- function(data, cl, alt){
	data1 <- data[,cl == 1]
	data2 <- data[,cl == 0]

	if(alt == 'greater'){
		rowMax <- apply(data1,1,max)
		pastThreshold <- data2 - rowMax
	}
	if(alt == 'less'){
		rowMin <- apply(data1,1,min)
		pastThreshold <- rowMin - data2
	}
	
	pastThreshold[pastThreshold <= 0] <- NA

	return(pastThreshold)
}


#Scoring function. Current implementation allows 'tanh' and 'quad' options. ModF specifies the parameter used within the scoring function.
cutf <- function(x, modF=2, modType='tanh'){
	if(modType == 'tanh'){
		return(tanh((x-1)*(1+modF)) + 1)
	}
	if(modType == 'quad'){
		return(x^modF)
	}
}


#Permutation based (analytical) p-value calculation
test.tabulation <- function(n1, n2){
    logProbaStat = lchoose(n = n1+n2-(1:(n2+1)), k=n1-1) - lchoose(n= n1+n2, k=n1);
    probaStat = exp(logProbaStat);
    pvalues = cumsum(probaStat[(n2+1):1])[(n2+1):1];
    return(pvalues);
}


#convenience function to get candidate features, ordered by their DIDS score and filtered by p-value
#didsscore: DIDS scores (i.e. resulting from call to DIDSscore)
#pvalcutiff: the p-value used to filter the results by (i.e. anything with a larger p-value than this will be filtered out)
getCandidates <- function(didsscores, pvalcutoff=0.05){
	candidateIndex <- order(abs(didsscores$score), decreasing=T)
	candidates <- data.frame(index=candidateIndex, score=didsscores$score[candidateIndex], pval=didsscores$pval[candidateIndex])
	candidates <- candidates[candidates$pval < pvalcutoff, ]
}



#convenience function to make a matrix of samples that are aberrantly expressed in the cases, compared to the controls
#data: matrix of data; genes per row, samples per column (i.e. data on which DIDS was run)
#cl: class vector; 1 = controls/responder/sensitive; 0 = cases/non-responder/resistant
#candidatesgenes: the genes to include in the return matrix
#alt: specify wether candidate genes are aberrantly over- or under-expressed (i.e. alt='greater' or alt='less' respectively)
explainedCases <- function(data, cl, candidategenes, alt='greater', discretize=T){
	dataCandidates <- data[candidategenes, ]
	pastThreshold <- pastThreshold(dataCandidates, cl, alt)

	if(discretize){
		pastThreshold[!is.na(pastThreshold)] <- 1
	}
	pastThreshold[is.na(pastThreshold)] <- 0

	colnames(pastThreshold) <- colnames(dataCandidates)[cl == 0]
	rownames(pastThreshold) <- rownames(dataCandidates)

	return(pastThreshold)
}
