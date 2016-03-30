#plot expression values for a particular gene, separating by class
#data: data in eSet format
#class: class labels
#sorted: sort up or down (towards class-separator-middle)
#cols: class colors
#colCode: color per sample 
plotGroups <- function(data, class, gene=NULL, densPlot=T, binclass=T, sorted=0, groupNames=c('group1', 'group2'), main=NULL, ylim=c(0,16), ylab='Expression', cols=NULL, colCode=NULL, minkw=NULL, selection=NULL, ...){
	
	
	if(!is.null(selection)){
		data <- data[,selection]
		class <- class[selection]
	}
	
	if(!is.null(gene)){
		if(!is.numeric(gene)){
			geneIndex <- grep(paste('^',gene,'$', sep=""), fData(data)$symbol)
			stopifnot(length(geneIndex) == 1)
			geneIndex <- geneIndex[1]
			geneName <- gene		
			#todo: allow multiple plots for all probes for gene	
		}else{
			geneIndex <- gene
			geneName <- fData(data)$symbol[gene]
		}
		
		data <- exprs(data)[geneIndex,,drop=F]
	}
	else{
		geneName = ''
	}

	#make sure class data is complete
	data <- data[,!is.na(class),drop=F]
	
	completeClass <- class[!is.na(class)]
	if((length(class) - length(completeClass)) > 0){
		print(paste("Removing",(length(class) - length(completeClass)),"incomplete cases"))
		class <- completeClass
	}
	
	if(binclass){
		index <- vector(mode="list", length=2)
		index[[1]] <- which(class)
		index[[2]] <- which(!class)
		totalNr <- ncol(data)
	}
	else{
		distinctClasses <- levels(as.factor(class))
		index <- vector(mode="list", length=length(distinctClasses))
		totalNr <- 0
		
		for(i in 1:length(distinctClasses)){
			index[[i]] <- which(class == distinctClasses[i])
			totalNr <- totalNr + length(index[[i]])
		}
	}
	
	if(is.null(colCode)){ colCode <- rep(NA, ncol(data)) }
	if(is.null(cols)){ cols <- 1:length(index) }
	
	if(is.null(main)){
		main=paste('Expression of',geneName,'for',paste(groupNames, collapse=' vs '))
	}
	
	if(densPlot){
		layout(matrix(c(1,1,1,1,1,1,1,2), ncol=4, byrow=T))
	}
	
	curX <- 0
	for(i in 1:length(index)){
		futX <- length(index[[i]]) + curX
		plotData <- data[,index[[i]]]
		
		if(((i+sorted) %% 2) > 0){
			sortOrder <- F
		}else{
			sortOrder <- T
		}
		
		if(sorted >= 0){
			sortedOrder <- order(plotData, decreasing=sortOrder)
			plotData <- plotData[sortedOrder]
		}		
		
		if(i == 1){
			plot((curX+1):futX, plotData, xlim=c(1,totalNr+2), ylim=ylim, col=cols[i], xlab='Samples', ylab=ylab, main=main, ...)
		}
		else{
			points((curX+1):futX, plotData, col=cols[i], ...)
		}
		
		if(sum(!is.na(colCode))>0){
			points((curX+1):futX, plotData, col=colCode[index[[i]]][sortedOrder], pch=20)
		}

		abline(v=curX, col='grey')
		
		points(curX, mean(data[,index[[i]]]), pch=20, col=i)
		text(curX,ylim[1], groupNames[i], pos=4)
	
		curX <- futX + 2
	}
	
	if(length(data[,index[[1]]]) > 1 & length(data[,index[[2]]]) > 1 & binclass==T){
		stat <- signif(wilcox.test(data[,index[[1]]], data[,index[[2]]])$p.value, 2)
		title(sub=paste('mann-whitney p-value:', stat))
	}
	
	if(densPlot){
		dens1 <- density(data[,index[[1]]], na.rm=T)
		if(!is.null(minkw)){
			dens1 <- density(data[,index[[1]]], bw=minkw, na.rm=T)
		}
		
		plot(dens1, ann=F, axes=F, xlim=c(min(data, na.rm=T) - 1,max(data, na.rm=T) + 1))
		lines(density(data[,index[[2]]], bw=dens1$bw, na.rm=T), col='red')
		title(sub=paste('Bandwidth:', dens1$bw))
	}
}

#convenience function to create a pdf with the top x candidate genes
#data: data in eSet format
#class: class labels
#candidates: features to plot
#sorted: sort up or down (towards class-separator-middle)
#colCode: color per sample 
pdftop <- function(data, class, candidates, filename, top=NULL, ylim=c(0,15), binclass=T, groupNames=NULL, selection=NULL, sorted=0, colCode=NULL){
	pdf(file=filename)
	
	if(is.null(top)){
		top <- length(candidates)
	}
	for(i in 1:top){
		plotGroups(data, class, gene=candidates[i], ylim=ylim, binclass=binclass, selection=selection, groupNames=groupNames, sorted=sorted, colCode=colCode)
	}
	dev.off()
}

#convenience function to plot explained cases (i.e. result from call to 'explainedCases' function)
#explainedCases: result from call to 'explainedCases' function
#featurenames: names for the features (i.e. rows)
#sampleNames: names for the samples (i.e. columns)
#additional parameters that can be parsed by the R 'plot' function can be passed to this function
plotExplainedCases <- function(explainedCases, featureNames=rownames(explainedCases), sampleNames=colnames(explainedCases), ...) {
	heatmap(explainedCases, scale='none', col=c('white', 'blue'), ...)
}