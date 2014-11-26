# DESeqTools.R

# DESeq -- differential expression analysis of digital gene expression data

#  	Simon Anders, Wolfgang Huber


DESeq.DiffExpress <- function( fnames, fids, groupSet, targetGroup=sort(groupSet)[1], geneColumn="GENE_ID", 
			intensityColumn="READS_M", keepIntergenics=FALSE, 
			minimumRPKM=1, missingGenes="fill", extraColumn=NULL,
			average.FUN=sqrtmean, wt.folds=1, wt.pvalues=1, wt.dists=1, 
			...) {

	# turn the set of transcript files into one matrix
	m <- expressionFileSetToMatrix( fnames=fnames, fids=fids, geneColumn=geneColumn,
			intensityColumn=intensityColumn, missingGenes=missingGenes)
	if ( ! is.null( extraColumn)) {
		mExtra <- expressionFileSetToMatrix( fnames=fnames, fids=fids, geneColumn=geneColumn,
				intensityColumn=extraColumn, missingGenes=missingGenes)
	}


	# drop non-genes...
	if ( ! keepIntergenics) {
		drops <- grep( "(ng)", rownames(m), fixed=T)
		if ( length(drops) > 0) {
			m <- m[ -drops, ]
			if ( ! is.null( extraColumn)) mExtra <- mExtra[ -drops, ]
		}
		nonGenes <- subset( getCurrentGeneMap(), REAL_G == FALSE)$GENE_ID
		drops <- which( rownames(m) %in% nonGenes)
		if ( length(drops) > 0) {
			m <- m[ -drops, ]
			if ( ! is.null( extraColumn)) mExtra <- mExtra[ -drops, ]
		}
	}

	# the column flags for DESeq are to end up as 1,2, where 1 is the one we want...
	cl <- rep( 2, times=length(groupSet))
	cl[ groupSet == targetGroup] <- 1
	myGrpNames <- rep( targetGroup, times=length(groupSet))
	myGrpNames[ cl == 2] <- notTargetGroup <- paste( "Not", targetGroup, sep=" ")
	if ( ! is.null( extraColumn)) {
		clExtra <- cl
	}

	# the dispersion has to be estimated differently if there aren't replicates
	specialDispersion <- FALSE
	if ( any( table(cl) < 2)) specialDispersion <- TRUE

	mUse <- matrix( as.integer(m), nrow=nrow(m), ncol=ncol(m))
	colnames(mUse) <- colnames(m)
	rownames(mUse) <- rownames(m)

	require( DESeq)
	NG <- nrow( m)

	# call DESeq,  it's a multi step process
	ans <- newCountDataSet( countData=mUse, conditions=cl)
	ans <- estimateSizeFactors( ans)
	if (specialDispersion) {
		ans <- estimateDispersions( ans, method="blind", sharingMode="fit-only", fitType="local")
	} else {
		ans <- estimateDispersions( ans, fitType="local")
	}
	deseqOut <- nbinomTest( ans, 2, 1)

	# DESeq throws away rows that were all zero, and may do divide by zero
	gnames <- rownames(mUse)
	deseqGenes <- deseqOut[ ,1]
	gPtr <- match( gnames, deseqGenes, nomatch=0)
	gprod <- gene2ProductAllSpecies( gnames)
	if ( any( gprod == "")) {
		tmp <- read.delim( fnames[1], as.is=T)
		gnams <- tmp[[ geneColumn]]
		gpros <- tmp[[ grep( "product", tolower(colnames(tmp)))]]
		needs <- which( gprod == "")
		where <- match( gnames[needs], gnams, nomatch=0)
		gprod[ needs[ where > 0]] <- gpros[ where]
	}

	v1 <- v2 <- foldOut <- rep.int( 0, length(gnames))
	v1[ gPtr > 0] <- deseqOut$baseMeanA[ gPtr]
	v2[ gPtr > 0] <- deseqOut$baseMeanB[ gPtr]
	foldOut <- log2( (v2+minimumRPKM) / (v1+minimumRPKM))
	pvalOut <- rep.int( 1, length(gnames))
	pvalOut[ gPtr > 0] <- deseqOut$pval[ gPtr]
	pvalOut[ is.na( pvalOut)] <- 1

	out <- data.frame( gnames, gprod, foldOut, pvalOut, v2, v1,
			stringsAsFactors=F)
	colnames(out) <- c( "GENE_ID", "PRODUCT", "LOG2FOLD", "PVALUE", 
			targetGroup, notTargetGroup)

	if ( ! is.null( extraColumn)) {
		avgExtra1 <- apply( mExtra[ , which( cl == 1)], MARGIN=1, FUN=average.FUN)
		avgExtra2 <- apply( mExtra[ , which( cl == 2)], MARGIN=1, FUN=average.FUN)
		out$AVG_EXTRA1 <- avgExtra1
		out$AVG_EXTRA2 <- avgExtra2
	}
	
	# the fold change data is biased by raw read counts, ignore it!!
	ord <- diffExpressRankOrder( out$LOG2FOLD, out$PVALUE, wt.folds, wt.pvalues)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)

	return( out)
}

