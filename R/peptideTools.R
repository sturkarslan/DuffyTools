# peptideTools.R -  peptide based tools


DNAtoBestPeptide <- function( dnaSet, clipAtStop=FALSE) {

	out <- sapply( as.character( dnaSet), function( dna) {
			peps <- DNAtoAA.fast( dna, clipAtStop=clipAtStop)
			if (clipAtStop) {
				return( peps[ which.max( nchar(peps))])
			} else {
				#terms <- unlist( strsplit( peps, split=STOP_CODON, fixed=TRUE), use.names=F)
				terms <- unlist( strsplit( peps, split="\\*|\\?", fixed=FALSE), use.names=F)
				return( terms[ which.max( nchar(terms))])
			}
		}, USE.NAMES=FALSE)

	return( out)
}


peptideChopper <- function( proteinSet, len=15, overlap=9, verbose=TRUE) {

	# find all the peptides of length "len" with overlap "overlap" in the
	# given protein sequences...

	if (verbose) cat( "\n\n-----------------------\nPeptide Length: ", len, 
				"\tOverlap: ", overlap)

	# visit each protein one at a time
	out <- list()
	outCounts <- vector()

	for (i in 1:length( proteinSet)) {
		pseq <- proteinSet[i]
		inow <- 1
		ilast <- nchar( pseq)
		if ( ilast < len) {
			warning( paste( "Protein smaller than peptide length:", gSet[i]))
			next
		}
		npeps <- 0
		pepSet <- pepStart <- pepEnd <- pepOrdinal <- vector()
		perfectEnding <- FALSE

		# step along the protein, accumulating peptides
		repeat {
			ito <- min( c( (inow+len-1), ilast))
			pepNow <- substr( pseq, inow, ito)
			nNow <- nchar( pepNow)
			# if its not the full peptide length, then we just stepped off the end
			if ( nNow < len) break
			# add this peptide to the growing list (no worry about duplicates yet)
			npeps <- npeps + 1
			pepSet[npeps] <- pepNow
			pepStart[npeps] <- inow
			pepEnd[npeps] <- ito
			pepOrdinal[npeps] <- npeps
			# watch for the special case of a peptide ending exactly at the end of the protein
			perfectEnding <- (ito == ilast)
			# step ahead, with the right overlap of the previous peptide
			inow <- ito - overlap + 1
		}
		# after the end, make one last peptide that has greater overlap, just to
		# make sure we end on the last a.a. of the protein
		if ( ! perfectEnding) {
			pepNow <- substr( pseq, (ilast-len+1), ilast)
			npeps <- npeps + 1
			pepSet[npeps] <- pepNow
			pepStart[npeps] <- (ilast-len+1)
			pepEnd[npeps] <- ilast
			pepOrdinal[npeps] <- npeps
		}
		# eliminate duplicates
		nBefore <- length( pepSet)
		duplicates <- which( duplicated( pepSet))
		if ( length( duplicates) > 0) {
			pepSet <- pepSet[ -duplicates]
			pepStart <- pepStart[ -duplicates]
			pepEnd <- pepEnd[ -duplicates]
			pepOrdinal <- pepOrdinal[ -duplicates]
		}
		nAfter <- length( pepSet)
		if (verbose) cat( "\n\tGene: ", names(proteinSet)[i], "\tN_Peptides: ", 
					nAfter,"\tN_Duplicates_Removed:", nBefore-nAfter)
		
		# make a small data frame that has "where info"
		myDF <- data.frame( pepSet, pepStart, pepEnd, pepOrdinal, stringsAsFactors=FALSE)
		colnames( myDF) <- c("Peptide", "Position", "End", "Ordinal")
		out[[i]] <- myDF
		names( out)[i] <- names(proteinSet)[i]
		outCounts[i] <- nAfter
	}
	return( list( "peptideSets"=out, "peptideCounts"=outCounts) )
}


peptideTable <- function( gSet, len=15, overlap=9, verbose=TRUE) {

	out <- peptideChopper( gSet, len, overlap, verbose=verbose)
	# just the list of data frames
	ans <- out$peptideSets

	# build a data.frame...
	nrows <- 0
	ncols <- length( ans)
	totalPeps <- 0
	for ( i in 1:ncols) {
		# "chopper" returns DFs now, not just vectors...
		myDF <- ans[[i]]
		nr <- length( myDF$Peptide)
		nrows <- max( c(nrows, nr))
		totalPeps <- totalPeps + nr
	}
	peps <- matrix( "", nrow=nrows, ncol=ncols)
	allPepNames <- vector()
	for ( i in 1:ncols) {
		myDF <- ans[[i]]
		pset <- myDF$Peptide
		np <- length( pset)
		peps[ 1:np, i] <- pset
		allPepNames <- union( allPepNames, pset)
	}
	colnames(peps) <- names( ans)
	write.table( peps, file="peptideTable.tsv", sep="\t", quote=FALSE)
	cat("\n\nSummary:  Length=", len, "\tOverlap=", overlap)
	cat("\n\tTotal Peptides=", totalPeps, "\tTotal Unique Peptides=", length(allPepNames),"\n")
}


peptideTableSet <- function( gSet, len=15, overlap=c(7:9)) {

	# call the peptide table and chopper iteratively to make a data frame of counts
	cnts <- matrix( nrow=length( gSet), ncol=length(overlap))

	for( i in 1:length(overlap)) {
		thislap <- overlap[i]
		out <- peptideChopper( gSet, len, thislap, FALSE)
		# we want the counts only
		cnts[ ,i] <- out$peptideCounts
	}

	myDF <- data.frame( gSet, cnts, stringsAsFactors=FALSE)
	colnames( myDF) <- c( "Gene", paste( "Overlap_", overlap, sep=""))
	rownames( myDF) <- 1:nrow( myDF)
	write.table( myDF, file="peptideOverview.txt", sep="\t", quote=FALSE)
}


peptideViewer <- function( gene, AAoffset=0, len=15, overlap=9, overlayDF=NULL, asPNG=FALSE) {

	# draw where the peptides land...
	# protein <<- extractProteinSequenceFromWeb( gene, verbose=FALSE)
	protein <<- gene2Protein( gene)
	out <<- peptideChopper( gene, len=len, lap=lap, force=TRUE)
	ans <- out$peptideSets

	if (asPNG) {
		jpeg( file=paste( gene, ".jpg", sep=""), width=1000, height=650, pointsize=12)
	}

	# plot setup
	myCex <- 1.0
	par( family="mono")
	nc <- nchar( protein[1])
	ncDraw <- par("pin") / par("cin") * 2

	# try to overlay other peptides to see where they land...  Christian's for now
	overlay <- FALSE
	if ( ! is.null( overlayDF)) {
		otherPeps <- subset( overlayDF, GeneID==gene)
		if ( nrow( otherPeps) > 0) {
			overlay <- TRUE
			otherStarts <- vector()
			for( i in 1:nrow( otherPeps)) otherStarts[i] <- regexpr( otherPeps$Sequence[i], protein, fixed=TRUE)
			missing <- which( otherStarts < 0)
			if (length( missing) > 0) cat( "\nN_overlay peptides not matching protein", length(missing))
		}
	}

	# total number of plots based on how many aa we can fit across X axis
	nPlots <- ceiling( nc / ncDraw)[1]
	nPlotsPerScreen <- 3
	par(mfcol=c(nPlotsPerScreen,1))
	for ( iplot in 1:nPlots) {
		thisXlim <- c( ((iplot-1)*ncDraw[1]), iplot*ncDraw[1]) + AAoffset
		plot( 0,0, xlim=thisXlim, ylim=c(0,5),type="n", xlab=paste( gene, "  amino acid position"), ylab="")
	
		# blindly draw the whole protein at the bottom
		text( (1:nc), rep(0,times=nc), labels=unlist(strsplit(protein,"")), adj=0, cex=myCex, col=1)
		xwid <- 1.0
		ygap <- 0.7

		# step through the peptides and draw them
		maxDeep <- floor( len / (len-lap)) + 1
		myDF <- ans[[1]]
		for ( i in 1:nrow(myDF)) {
			pep <- myDF$Peptide[i]
			ifrom <- myDF$Position[i]
			ito <- myDF$End[i]
			ordinal <- myDF$Ordinal[i]
			curDeep <- ( (ordinal-1) %% maxDeep) + 1
			text( (((ifrom-1)*xwid)+1), (curDeep*ygap), labels=pep, adj=0, cex=myCex, col=4)
			if ( ordinal %% 10 == 0) {
				text( (((ifrom-1)*xwid)+2), ((curDeep+0.5)*ygap), labels=ordinal, adj=0, cex=myCex*0.95, col=2)
			}
		}

		# now try to overlay the 'others'
		if ( overlay) {
			overSet <- which( (otherStarts > 0) & (otherStarts >= (thisXlim[1]-len)) & (otherStarts <= thisXlim[2]))
			if ( length( overSet) > 0) {
				text( (((otherStarts[overSet]-1)*xwid)+1), ((maxDeep+1)*ygap), 
					labels=otherPeps$Sequence[ overSet], adj=0, cex=myCex, col=2)
				text( (((otherStarts[overSet]-1)*xwid)+2), ((maxDeep+2)*ygap), 
					labels=otherPeps$PeptideID[overSet], adj=0, cex=myCex, col=2)
			}
		}

		# do we need to pause for screen capture?
		if ( asPNG) {

		} else {
			if ( iplot %% nPlotsPerScreen == 0 && iplot < nPlots) {
				cat( "\npausing for screeen capture.  Click on plot to continue.")
				locator(1)
			}
		}
	}
	if ( asPNG) dev.off()
}


peptideRedundancySniffer <- function() {

	# this will see if we can save many peptides from being ordered, by checking
	# for big partial match from one end of peptides
	cat( "\nReading in peptide table:     peptideTable.txt")
	pepTable <- read.delim( "peptideTable.txt", as.is=TRUE)

	# gather up all peptides as one big set
	pepSet <- vector()
	for ( i in 1:ncol(pepTable)) {
		oneSet <- pepTable[ ,i]
		# because the table is ragged, some are empty slots
		good <- which( oneSet != "")
		cat( "\nReading Gene: ", colnames(pepTable)[i], "\tN_peptides in: ", length( good))
		pepSet <- append( pepSet, oneSet[ good])
	}
	cat( "\nTotal Peptides read in: ", length( pepSet))
	cat( "\nTotal Unique Peptides (full length): ", length( unique( pepSet)))

	fullLen <- nchar( pepSet[1])
	for ( iLess in 1:3) {
		iend <- fullLen - iLess
		smlSet <- substr( pepSet, 1, iend)
		cat( "\n\nUnique over N-terminal ",iend, "-mer:  ", length( unique( smlSet)), sep="")
	}
	for ( iLess in 1:3) {
		iend <- fullLen
		ibeg <- 1 + iLess
		smlSet <- substr( pepSet, ibeg, iend)
		cat( "\n\nUnique over C-terminal ",(fullLen-iLess), "-mer:  ", length( unique( smlSet)), sep="")
	}
}


peptide.mergeOverlaps <- function( peptides, max.tail=5, min.count=2) {

	if ( is.data.frame( peptides)) {
		if ( ! all( colnames( peptides) == c( "Peptide", "Count"))) 
				stop( "Expected data frame with {Peptide, Count} columns")
		ord <- order( peptides$Count, decreasing=TRUE)
		allPeps <- peptides$Peptide[ord]
		allCnts <- peptides$Count[ord]
	} else {
		allPeps <- peptides
		allCnts <- rep.int( 1, length(peptides))
	}

	cat( "\nGiven Peptides:  ", N <- length( allPeps))

	tapply( 1:N, factor( allPeps), FUN=function(x) {
			if ( length(x) < 2) return()
			newcnt <- sum( allCnts[x])
			allCnts[ x[1]] <<- newcnt
			allCnts[ x[2:length(x)]] <<- 0
			return()
		}, simplify=FALSE)

	peps <- allPeps[ allCnts > 0]
	cnts <- allCnts[ allCnts > 0]
	cat( "\nUnique Peptides: ", length( peps))

	# only keep the ones with at least K reads to support it
	keep <- which( cnts >= min.count)
	peps <- peps[ keep]
	cnts <- cnts[ keep]
	len <- nchar( peps)
	N <- length( peps)
	cat( "\nCounts >= ",min.count,":  ", N)


	# build a set of A-Z 'first AA' sets
	first <- substr( peps, 1,1)
	letterSets <- lapply( LETTERS, function(x) which( first == x))
	names( letterSets) <- LETTERS

	cat( "\n")
	savedIters <<- vector( mode="list")

	for ( iter in 1:max.tail) {
		nmerge <- 0
		for ( i in 1:N) {

			if ( len[i] < iter) next
			# extract the suffix from one peptide, and see all it hits
			mykey <- substr( peps[i], 1, 1)
			prefix <- substr( peps[i], 1, iter)
			suffix <- substr( peps[i], iter+1, len[i])
			theirkey <- substr( suffix, 1,1)
			totest <- letterSets[[ theirkey]]
			hits <- grep( suffix, peps[totest], fixed=T)
			#targets <- substr( peps[totest], 1, nchar(suffix))
			#hits <- which( targets == suffix)
			if ( length(hits) < 1) next
			hits <- totest[ hits]
			hits <- setdiff( hits, i)
			if ( length(hits) < 1) next

			# for each other peptide that starts with my suffix, prepend my prefix
			for ( j in hits) {
				loc <- regexpr( suffix, peps[j], fixed=TRUE)
				if ( loc != 1) next
				newpep <- paste( prefix, peps[j], sep="")
				letterSets[[ theirkey]] <- setdiff( letterSets[[ theirkey]], j)
				letterSets[[ mykey]] <- c( letterSets[[ mykey]], j)
				peps[j] <- newpep
				len[j] <- nchar(newpep)
				cnts[j] <- cnts[j] + cnts[i]
				nmerge <- nmerge + 1
			}

			# remove me fron the set
			letterSets[[mykey]] <- setdiff( letterSets[[mykey]], i)
			peps[i] <- ""
			len[i] <- 0
			cnts[i] <- 0
			if ( i %% 100 == 0) cat( "\r", i, nmerge)
		}
		cat( "\nIteration: ", iter, "\nTable of peptide lengths:\n")
		print( ans <- table( nchar( peps)))
		cat( "\nTotal:  ", sum( ans), "\n")
		savedIters[[iter]] <<- peps
	}

	keep <- which( len > 0 & cnts > 0)
	out <- data.frame( "Peptide"=peps[keep], "Count"=cnts[keep], stringsAsFactors=FALSE)
	ord <- order( out$Peptide, -out$Count)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)

	return( out)
}
