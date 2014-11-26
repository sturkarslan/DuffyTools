# metaRank.R -- combine several ranked files of genes

 
metaRanks <- function( fnames, fids, weightset=rep(1, length(fnames)), 
			geneColumn="GENE_ID", valueColumn="LOG2FOLD", 
			pvalueColumn="PVALUE", productColumn="PRODUCT", sep="\t",
			rank.average.FUN=sqrtmean, value.average.FUN=mean,
			keepIntergenics=FALSE, missingGenes=c("drop", "fill", "na"), 
			missingValue=0, naDropPercent=0.5, nFDRsimulations=0) {

	missingGenes <- match.arg( missingGenes)

	cat( "\nReading in files...")
	allDF <- vector( "mode"="list")
	allGenes <- vector()
	nfiles <- 0
	missing <- vector()
	for (i in 1:length(fnames)) {
		filename <- fnames[i]
		if ( ! file.exists( filename)) {
			cat( "\nNot found: ", filename)
			missing <- c( missing, i)
			next
		}
		tmp <- read.delim( filename, as.is=T, sep=sep)
		if ( ! geneColumn %in% colnames(tmp)) {
			cat( "\nGeneID column not found.   File=", filename)
			next
		}
		nfiles <- nfiles + 1
		if ( ! keepIntergenics) {
			isNG <- grep( "(ng)", tmp[[ geneColumn]], fixed=T)
			if ( length(isNG) > 0) tmp <- tmp[ -isNG, ]
			nonGenes <- subset( getCurrentGeneMap(), REAL_G == FALSE)$GENE_ID
			drops <- which( tmp[[ geneColumn]] %in% nonGenes)
			if ( length(drops) > 0) tmp <- tmp[ -isNG, ]
		}
		allDF[[ nfiles]] <- tmp
		allGenes <- base::union( allGenes, tmp[[geneColumn]])
		cat( "\n",filename, "\tN_Genes: ", nrow(tmp))
	}

	# watch for Human gene name complications...
	#trapHuman <- FALSE
	#if ( length( grep( ".Hs.", fnames, fixed=T)) > 0) {
	#	allGenes <- sub( ":.+", "", allGenes)
	#	trapHuman <- TRUE
	#}
	allG <- sort( unique( allGenes))
	NG <- length( allG)

	if ( length( missing) > 0) {
		fnames <- fnames[ -missing]
		fids <- fids[ -missing]
	}
	cat( "\nN_Files: ", nfiles, "   N_Genes: ", NG)
	
	# get the rank order in all datasets
	rankM <- foldM <- pvalM <- matrix( NA, nrow=NG, ncol=nfiles)
	rownames(rankM) <- rownames(foldM) <- rownames(pvalM) <- allG
	colnames(rankM) <- colnames(foldM) <- colnames(pvalM) <- fids
	allProds <- rep( "", NG)

	for( i in 1:nfiles) {
		thisDF <- allDF[[i]]
		theseGenes <- thisDF[[geneColumn]]
		#if (trapHuman) theseGenes <- sub( ":.+", "", theseGenes)
		where <- match( allG, theseGenes, nomatch=0)
		rankM[ where > 0, i] <- where[ where > 0]
		logFoldColumn <- grep( valueColumn, colnames(thisDF))
		if ( length( logFoldColumn) > 0) {
			foldM[ where > 0, i] <- thisDF[[ logFoldColumn[1]]][where]
		}
		pvalColumn <- grep( pvalueColumn, colnames(thisDF))
		if ( length( pvalColumn) > 0) {
			pvalM[ where > 0, i] <- thisDF[[ pvalColumn[1]]][where]
		}
		hasPROD <- ( productColumn %in% colnames(thisDF))
		if ( hasPROD) {
			allProds[ where > 0] <- ifelse( allProds[where > 0] == "", 
							thisDF[[ productColumn]][where], allProds[ where > 0])
		}
	}

	# check for missing genes
	whoNA <- vector()
	for ( i in 1:nfiles) {
		who <- which( is.na( rankM[ ,i]))
		if ( length(who) > 0) {
			whoNA <- sort( unique( c( whoNA, who)))
			if ( missingGenes == "fill") {
				rankM[ who, i] <- nrow(rankM)
				foldM[ who, i] <- missingValue
				pvalM[ who, i] <- 1.0
			}
		}
	}
	if ( missingGenes == "drop" && length(whoNA) > 0) {
		rankM <- rankM[ -whoNA, ]
		foldM <- foldM[ -whoNA, ]
		pvalM <- pvalM[ -whoNA, ]
		allG <- allG[ -whoNA]
		allProds <- allProds[ -whoNA]
		NG <- nrow(rankM)
	}

	# if 'na', and a gene row has too many na's, drop the whole thing
	if ( missingGenes == "na") {
		nna <- apply( rankM, 1, function(x) sum(is.na(x)))
		whoNA <- which( nna > (nfiles * naDropPercent))
		if ( length( whoNA) > 0) {
			rankM <- rankM[ -whoNA, ]
			foldM <- foldM[ -whoNA, ]
			pvalM <- pvalM[ -whoNA, ]
			allG <- allG[ -whoNA]
			allProds <- allProds[ -whoNA]
			NG <- nrow(rankM)
		}
	}
	cat( "    after 'drop' and/or 'na' filtering: ", NG)

	cat( "\nAveraging Ranks...")
	if ( all( weightset == 1)) {
		avgRank <- apply( rankM, MARGIN=1, FUN=rank.average.FUN, na.rm=T)
	} else {
		if ( identical( rank.average.FUN, logmean)) {
			avgRank <- apply( rankM, MARGIN=1, FUN=function(x) {
						numerator <- sum( log2(x) * weightset)
						denom <- sum( weightset)
						return( 2 ^ (numerator / denom))
					})
		} else if ( identical( rank.average.FUN, sqrtmean)) {
			avgRank <- apply( rankM, MARGIN=1, FUN=function(x) {
						numerator <- sum( sqrt(x) * weightset)
						denom <- sum( weightset)
						return( (numerator / denom) ^ 2)
					})
		} else {
			avgRank <- apply( rankM, MARGIN=1, FUN=function(x) {
						numerator <- sum( x * weightset)
						denom <- sum( weightset)
						return( numerator / denom)
					})
		}
	}

	avgFold <- apply( foldM, MARGIN=1, FUN=value.average.FUN, na.rm=T)
	avgPval <- apply( pvalM, MARGIN=1, FUN=logmean, na.rm=T)
	out <- data.frame( allG, allProds, avgFold, avgPval, avgRank, rankM, 
			stringsAsFactors=FALSE)
	colnames(out) <- c( "GENE_ID", "PRODUCT", valueColumn, "AVG_PVALUE", "AVG_RANK", colnames( rankM))

	# do the final ordering by Average Rank
	ord <- order( out$AVG_RANK, out$AVG_PVALUE, -out[[ valueColumn]])
	out <- out[ ord, ]
	rownames( out) <- 1:nrow(out)

	# do a simulation of random permutations of these ranks
	if ( nFDRsimulations > 0) {
		simM <- rankM
		randomAvgRank <- vector()
		cat( "  estimating FDR..")
		for ( k in 1:nSimulations) {
			for ( i in 1:nfiles) simM[ , i] <- sample( simM[ ,i])
			randomAvgRank <- c( randomAvgRank, 
						apply( simM, MARGIN=1, FUN=rank.average.FUN, na.rm=T))
		}
		# with this pool of 'by-chance average ranks, we can estimate the likelihood of ours
		randomAvgRank <- sort( randomAvgRank)
		avgRank <- out$AVG_RANK
		myLocs <- findInterval( avgRank * 1.00001, randomAvgRank)
		myLocs <- ifelse( myLocs > 0, myLocs - 1, 0)
		myEvalue <- myLocs / nSimulations
		myFPrate <- myEvalue / (1:nrow(out))
		myFPrate <- ifelse( myFPrate > 1, 1, myFPrate)
		out$E_VALUE <- myEvalue
		out$FP_RATE <- myFPrate
	}

	# correlation test...
	ccM <- matrix( NA, nrow=nfiles, ncol=nfiles)
	for( i in 1:(nfiles-1)) {
	for( j in (i+1):nfiles) {
		thisCC <- cor( rankM[ ,i], rankM[ ,j], use="complete")
		ccM[i,j] <- ccM[j,i] <- thisCC
	}}
	colnames( ccM) <- colnames(rankM)
	cc <- apply( ccM, MARGIN=2, mean, na.rm=T)
	metaRankCC <<- cc

	cat( "\nRank Correlations:\n")
	print( sort( cc, decreasing=T))

	return( out)
}


metaRank2html <- function( tbl, fileout="metaRanks.html", title="", maxRows=100, 
			valueColumn="LOG2FOLD", ...) {

	# clean up any formatting...
	tbl[[ valueColumn]] <- formatC( tbl[[ valueColumn]], format="f", digits=3)
	tbl$AVG_PVALUE <- formatC( tbl$AVG_PVALUE, format="e", digits=3)
	tbl$AVG_RANK <- formatC( tbl$AVG_RANK, format="f", digits=2)
	colnames(tbl)[3:5] <- c( gsub( "_", " ", valueColumn), "Avg Pvalue", "Avg Rank")

	title <- paste( "Meta Ranks:  &nbsp; ", title)
	table2html( tbl, fileout=fileout, title=title, maxRows=maxRows, ...)
	return()
}


metaRank.data.frames <- function( df.list, weightset=rep(1, length(df.list)), 
			geneColumn="GENE_ID", valueColumn="LOG2FOLD", 
			pvalueColumn="PVALUE", productColumn="PRODUCT",
			rank.average.FUN=sqrtmean, value.average.FUN=mean,
			missingGenes=c("drop", "fill", "na"), missingValue=0,
			naDropPercent=0.5) {

	missingGenes <- match.arg( missingGenes)

	allDF <- vector( "mode"="list")
	nDF <- 0
	allGenes <- vector()
	for (i in 1:length(df.list)) {
		tmp <- df.list[[ i]]
		if ( is.null(tmp)) {
			cat( "\nData frame is NULL  DataFrame: ", i, "\t", names(df.list)[i])
			next
		}
		if ( nrow(tmp) < 1) {
			cat( "\nData frame is empty...  DataFrame: ", i, "\t", names(df.list)[i])
			next
		}
		if ( ! geneColumn %in% colnames(tmp)) {
			cat( "\nGeneID column not found.   DataFrame: ", i, "\t", names(df.list)[i])
			next
		}
		nDF <- nDF + 1
		allDF[[ nDF]] <- tmp
		names(allDF)[nDF] <- names(df.list)[i]
		allGenes <- base::union( allGenes, tmp[[geneColumn]])
		cat( "\n",i, "\tN_Genes: ", nrow(tmp))
	}
	fids <- names( allDF)

	allG <- sort( unique( allGenes))
	NG <- length( allG)
	cat( "\nN_DataFrames: ", nDF, "   N_Genes: ", NG)
	
	# get the rank order in all datasets
	rankM <- foldM <- pvalM <- matrix( NA, nrow=NG, ncol=nDF)
	rownames(rankM) <- rownames(foldM) <- rownames(pvalM) <- allG
	colnames(rankM) <- colnames(foldM) <- colnames(pvalM) <- fids
	allProds <- rep( "", NG)

	for( i in 1:nDF) {
		thisDF <- allDF[[i]]
		hasPROD <- ( productColumn %in% colnames(thisDF))
		theseGenes <- thisDF[[geneColumn]]
		where <- match( allG, theseGenes, nomatch=0)
		rankM[ where > 0, i] <- where[ where > 0]
		logFoldColumn <- grep( valueColumn, colnames(thisDF))
		if ( length( logFoldColumn) > 0) {
			foldM[ where > 0, i] <- thisDF[[ logFoldColumn[1]]][where]
		}
		pvalColumn <- grep( pvalueColumn, colnames(thisDF))
		if ( length( pvalColumn) > 0) {
			pvalM[ where > 0, i] <- thisDF[[ pvalColumn[1]]][where]
		}
		if (hasPROD) {
			allProds[ where > 0] <- ifelse( allProds[where > 0] == "", 
				thisDF[[ productColumn]][where], allProds[ where > 0])
		}
	}

	# check for missing genes
	whoNA <- vector()
	for ( i in 1:nDF) {
		who <- which( is.na( rankM[ ,i]))
		if ( length(who) > 0) {
			whoNA <- sort( unique( c( whoNA, who)))
			if ( missingGenes == "fill") {
				rankM[ who, i] <- nrow(rankM)
				foldM[ who, i] <- missingValue
				pvalM[ who, i] <- 1.0
			}
		}
	}
	if ( missingGenes == "drop" && length(whoNA) > 0) {
		rankM <- rankM[ -whoNA, ]
		foldM <- foldM[ -whoNA, ]
		pvalM <- pvalM[ -whoNA, ]
		allG <- allG[ -whoNA]
		allProds <- allProds[ -whoNA]
		NG <- nrow(rankM)
	}

	# if 'na', and a gene row has too many na's, drop the whole thing
	if ( missingGenes == "na") {
		nna <- apply( rankM, 1, function(x) sum(is.na(x)))
		whoNA <- which( nna > (nDF * naDropPercent))
		if ( length( whoNA) > 0) {
			rankM <- rankM[ -whoNA, ]
			foldM <- foldM[ -whoNA, ]
			pvalM <- pvalM[ -whoNA, ]
			allG <- allG[ -whoNA]
			allProds <- allProds[ -whoNA]
			NG <- nrow(rankM)
		}
	}
	cat( "    after 'drop' and/or 'na' filtering: ", NG)

	cat( "\nAveraging Ranks...")
	if ( all( weightset == 1)) {
		avgRank <- apply( rankM, MARGIN=1, FUN=rank.average.FUN, na.rm=T)
	} else {
		if ( identical( average.FUN, logmean)) {
			avgRank <- apply( rankM, MARGIN=1, FUN=function(x) {
						numerator <- sum( log2(x) * weightset)
						denom <- sum( weightset)
						return( 2 ^ (numerator / denom))
					})
		} else if ( identical( average.FUN, sqrtmean)) {
			avgRank <- apply( rankM, MARGIN=1, FUN=function(x) {
						numerator <- sum( sqrt(x) * weightset)
						denom <- sum( weightset)
						return( (numerator / denom) ^ 2)
					})
		} else {
			avgRank <- apply( rankM, MARGIN=1, FUN=function(x) {
						numerator <- sum( x * weightset)
						denom <- sum( weightset)
						return( numerator / denom)
					})
		}
	}

	avgFold <- apply( foldM, MARGIN=1, FUN=value.average.FUN, na.rm=T)
	avgPval <- apply( pvalM, MARGIN=1, FUN=logmean, na.rm=T)
	out <- data.frame( allG, allProds, avgFold, avgRank, avgPval, rankM, 
			stringsAsFactors=FALSE)
	colnames(out) <- c( "GENE_ID", "PRODUCT", valueColumn, "AVG_RANK", "AVG_PVALUE", colnames( rankM))

	# do the final ordering by Average Rank
	ord <- order( out$AVG_RANK, out$AVG_PVALUE, -out[[ valueColumn]])
	out <- out[ ord, ]
	rownames( out) <- 1:nrow(out)

	# correlation test...
	ccM <- matrix( NA, nrow=nDF, ncol=nDF)
	if ( nDF > 1) {
		for( i in 1:(nDF-1)) {
		for( j in (i+1):nDF) {
			thisCC <- cor( rankM[ ,i], rankM[ ,j], use="complete")
			ccM[i,j] <- ccM[j,i] <- thisCC
		}}
	}
	colnames( ccM) <- colnames(rankM)
	cc <- apply( ccM, MARGIN=2, mean, na.rm=T)
	metaRankCC <<- cc

	cat( "\nRank Correlations:\n")
	print( sort( cc, decreasing=T))

	return( out)
}
