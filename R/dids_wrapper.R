

didsScore <- function(eset, modF=2, modType='tanh', alt='two.sided') {
    stopifnot('condition' %in% colnames(phenoData(eset)))
    cl <- as.numeric(eset$condition == 'sensitive')
    didsScoreFromMatrix(data=exprs(eset), cl=cl, modF=modF, modType=modType, alt=alt)
}


didsScoreFromMatrix <- function(data, cl, modF=2, modType='tanh', alt='two.sided') {
    stopifnot(class(cl) == 'numeric')
    DIDSscore(data=data, cl=cl, modF=modF, modType=modType, alt=alt)
}


didsPlot <- function(eset, gene=NULL, geneSymbolCol='external_gene_id', 
                     groups=NULL, groupNames=c('Sensitive', 'Resistant'), 
                     densPlot=TRUE, binclass=TRUE, sorted=0, main=NULL, ylim=c(0,16), 
                     ylab='Expression', cols=NULL, colCode=NULL, minkw=NULL, selection=NULL, ...) {

    ## If group labels are provided, only check that these are logical so that
    ## they can be passed on to the plotGroups function. If no labels are provided
    ## we assume they can be derived from the condition column of the eset phenoData,
    ## which is assumed to contain 'sensitive' and 'resistant' values.
    if (!is.null(groups)) {
        stopifnot(class(groups) == 'logical')
    } else {
        stopifnot('condition' %in% colnames(phenoData(eset)))
        cl <- eset$condition == 'sensitive'
    }

    ## Inject 'symbol' column from specified column to conform with the
    ## expectancy of plotGroups that fData contains such a column containing
    ## the gene symbol of every gene in the dataset. If no geneSymbolCol is
    ## provided, we only check if the eset contains this column already 
    ## (and stop otherwise).
    if (geneSymbolCol != 'symbol') {
       stopifnot(geneSymbolCol %in% colnames(fData(eset)))
       fData(eset)$symbol <- fData(eset)[ , geneSymbolCol]
    } else {
        stopifnot('symbol' %in% colnames(fData(eset)))
    }

    ## Actually plot!
    plotGroups(data=eset, class=cl, gene=gene, densPlot=densPlot, binclass=binclass, sorted=sorted, 
               groupNames=groupNames, main=main, ylim=ylim, ylab=ylab, cols=cols, colCode=colCode, 
               minkw=minkw, selection=selection, ...)
}


expressionSetFromGCF <- function(countsPath, normalized=TRUE, log2=TRUE, samples=NULL,
                                 featCols=NULL, featIdCol='ensembl_gene_id') {
    ## Provide default GCF metadata columns if none provided
    if (is.null(featCols)) {
        featCols <- c('ensembl_gene_id', 'gene_biotype', 'chromosome_name',
                      'start_position', 'end_position', 'external_gene_id', 'description')
    }

    ## Load data from tab-delimited counts file
    data <- read.delim(countsPath, as.is=T, check.names=F)

    ## Extract feature data
    features <- data[ , featCols]
    rownames(features) <- features[ , featIdCol]
    featureData <- new('AnnotatedDataFrame', features)

    ## Extract count matrix
    counts <- as.matrix(data[ , !(colnames(data) %in% featCols)])
    rownames(counts) <- features[ , featIdCol]

    ## Normalize if required. Note that log2 implies normalization.
    ## DESeq2 is used for normalization and should be installed.
    if (normalized | log2) {
        counts <- normalizeCounts(counts)
        if (log2) counts <- log2(counts + 1)
    }

    ## Prepare sample data frame (if provided). Note that the count
    ## frame is subsetted to only include samples in the sample frame.
    if (!is.null(samples)) {
        phenoData <- new('AnnotatedDataFrame', samples)
        counts <- counts[ , rownames(samples)]
    } else {
        phenoData <-  annotatedDataFrameFrom(counts, byrow=FALSE)
    }

    ExpressionSet(assayData=counts, phenoData=phenoData, featureData=featureData)
}

normalizeCounts <- function(counts) {
  library(DESeq2)
  size.factors <- estimateSizeFactorsForMatrix(counts)
  t(t(counts) / size.factors)
}
