# samTools.R

# EdgeR -- differential expression analysis of digital gene expression data

#  	Mark Robinson, Davis McCarthy, Yunshun Chen, Gordon Smyth


EdgeR.DiffExpress <- function( fnames, fids, groupSet, targetGroup=sort(groupSet)[1], geneColumn="GENE_ID", 
			intensityColumn="READS_M", keepIntergenics=FALSE, 
			missingGenes="fill", extraColumn=NULL,
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

	# the column flags for EdgeR are to end up as 1,2, where 1 is the one we want...
	cl <- rep( 2, times=length(groupSet))
	cl[ groupSet == targetGroup] <- 1
	myGrpNames <- rep( targetGroup, times=length(groupSet))
	myGrpNames[ cl == 2] <- paste( "Not", targetGroup, sep=" ")
	if ( ! is.null( extraColumn)) {
		clExtra <- cl
	}

	# EdgeR will choke if less than two samples per group
	canDoDispersion <- TRUE
	dispersion <- "auto"
	if( sum( cl == 1) < 2) {
		canDoDispersion <- FALSE
	}
	if( sum( cl == 2) < 2) {
		canDoDispersion <- FALSE
	}

	grpFac <- factor( myGrpNames)
	grpNames <- levels(grpFac)

	mUse <- m

	require( edgeR)
	NG <- nrow( m)

	# call EdgeR,  it's a multi step process
	ans <- DGEList( counts=mUse, group=cl)
	ans <- calcNormFactors( ans)
	if ( canDoDispersion) {
		ans <- estimateCommonDisp( ans)
		ans <- estimateTagwiseDisp( ans)
	} else {
		dispersion <- 0.05
	}
	et <- exactTest( ans, pair=2:1, dispersion=dispersion)
	edgeRout <- topTags( et, n=NG)$table

	gnames <- rownames( edgeRout)
	where <- match( gnames, rownames(m), nomatch=0)
	mTmp <- m[ where, ]

	# calc the average for each group, and put it into 'this group' order
	avgM <- t( apply( mTmp, MARGIN=1, function(x) tapply(x, grpFac, FUN=average.FUN)))
	if ( colnames(avgM)[1] != targetGroup) avgM <- avgM[ ,c(2,1)]

	gprod <- gene2ProductAllSpecies( gnames)
	if ( any( gprod == "")) {
		tmp <- read.delim( fnames[1], as.is=T)
		gnams <- tmp[[ geneColumn]]
		gpros <- tmp[[ grep( "product", tolower(colnames(tmp)))]]
		needs <- which( gprod == "")
		where <- match( gnames[needs], gnams, nomatch=0)
		gprod[ needs[ where > 0]] <- gpros[ where]
	}

	out <- data.frame( gnames, gprod, edgeRout[ , c(1,3,4)], avgM, 
			stringsAsFactors=F)
	colnames(out) <- c( "GENE_ID", "PRODUCT", "LOG2FOLD", "PVALUE", 
			"FDR", colnames(avgM))

	if ( ! is.null( extraColumn)) {
		avgExtra1 <- apply( mExtra[ , which( cl == 1)], MARGIN=1, FUN=average.FUN)
		avgExtra2 <- apply( mExtra[ , which( cl == 2)], MARGIN=1, FUN=average.FUN)
		out$AVG_EXTRA1 <- avgExtra1
		out$AVG_EXTRA2 <- avgExtra2
	}
	
	# the fold change data is biased by raw read counts, ignore it!!
	ord <- diffExpressRankOrder( out$LOG2FOLD, out$PVALUE, wt.folds=0, wt.pvalues)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)

	return( out)
}

