# expressionFileSetToMatrix.R


`expressionFileSetToMatrix` <- function( fnames, fids, geneColumn=c( "GENE_ID", "GeneID"), 
					intensityColumn=c("INTENSITY","RPKM_M", "RANK"),
					missingGenes=c("na", "drop", "fill"), sep="\t",
					keepIntergenics=FALSE, verbose=FALSE) {

	missingGenes <- match.arg( missingGenes)

	# make sure we can read those files
	filesOK <- file.exists( fnames)
	if ( !all( filesOK)) {
		cat( "\nSome expression files not found:\n")
		print( fnames[ !filesOK])
		return(NULL)
	}
	nFiles <- length( fnames)

	# load each file in turn
	allData <- vector( mode="list")
	for( i in 1:nFiles) {
		tmp <- read.delim( fnames[i], as.is=T, sep=sep)
		if ( all( intensityColumn == "RANK") && (! ("RANK" %in% colnames(tmp)))) {
			tmp$RANK <- 1:nrow(tmp)
		}
		intenC <- base::match( intensityColumn, colnames(tmp), nomatch=0)
		if ( any( intenC > 0)) intenC <- intenC[ intenC > 0][1]
		geneC <- base::match( geneColumn, colnames(tmp), nomatch=0)
		if ( any( geneC > 0)) geneC <- geneC[ geneC > 0][1]
		if ( any( c(intenC,geneC) == 0)) {
			cat( "\nSome needed columns not found:   file: ", fnames[i],
				"\n  Expected: ", geneColumn, intensityColumn,
				"\n  Found:    ", colnames(tmp))
			return(NULL)
		}
		thisGenes <- tmp[[ geneC]]
		thisInten <- tmp[[ intenC]]
		smallDF <- data.frame( "GENE_ID"=thisGenes, "INTENSITY"=thisInten, stringsAsFactors=F)
		if ( ! keepIntergenics) {
			drops <- grep ( "(ng)", thisGenes, fixed=T)
			if ( length(drops) > 0) smallDF <- smallDF[ -drops, ]
			nonGenes <- subset( getCurrentGeneMap(), REAL_G == FALSE)$GENE_ID
			drops <- which( smallDF$GENE_ID %in% nonGenes)
			if ( length(drops) > 0) smallDF <- smallDF[ -drops, ]
		}
		allData[[i]] <- smallDF
		if ( i == 1) {
			allGenes <- smallDF$GENE_ID
		} else {
			if ( missingGenes == "drop") {
				nWas <- length(allGenes)
				allGenes <- base::intersect( allGenes, smallDF$GENE_ID)
				nNow <- length(allGenes)
				if ( nNow < nWas) cat( "\nSome geneIDs missing: ", basename(fnames[i]))
			} else {
				allGenes <- base::union( allGenes, smallDF$GENE_ID)
			}
		}
		if ( verbose) cat( "\nFile: ", i, basename(fnames[i]), "\tN_Genes: ", nrow(smallDF))
	}
	if (verbose) cat( "\n")

	allGenes <- base::sort( setdiff( allGenes, ""))
	nGenes <- length( allGenes)
	m <- matrix( NA, nrow=nGenes, ncol=nFiles)
	colnames(m) <- fids
	rownames(m) <- allGenes

	# now fill in the matrix
	smallestV <- max( thisInten, na.rm=T)
	for( i in 1:nFiles) {
		v <- rep( NA, times=nGenes)
		smallDF <- allData[[i]]
		thisGenes <- smallDF$GENE_ID
		thisInten <- smallDF$INTENSITY
		where <- base::match( allGenes, thisGenes, nomatch=0)
		v[ where > 0] <- thisInten[ where]
		m[ , i] <- v
		if ( ! any( intensityColumn == "RANK")) smallestV <- min( smallestV, thisInten, na.rm=T)
	}

	# if we kept all genes, fill in the holes with a small value
	if (missingGenes == "fill") {
		m[ is.na(m)] <- smallestV
	}

	return( m)
}
