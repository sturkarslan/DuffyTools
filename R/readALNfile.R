# readALNfile.R -- read up the result consensus alignment filr from Clustalw2 


readALNfile <- function( file, verbose=TRUE) {

	txt <- readLines( file)
	start <- 2
	end <- length(txt)
	gapLines <- which( txt == "")
	gapLines <- gapLines[ gapLines > start]
	starLines <- grep( "^    ", txt)

	# we can build a giant flat table of the alignments, finding group IDs as we go
	curLine <- start
	gptr <- 1
	NG <- 1
	big <- rep( "", times=NG)
	gnames <- rep( "", times=NG)
	starText <- ""
	starHeadPtr <- vector()

	while ( curLine <= end) {
		thisLine <- txt[ curLine]
		if ( curLine %in% starLines) {
			if ( length( starHeadPtr) > 0) {
				starStart <- median( starHeadPtr)
				thisStars <- substr( thisLine, starStart, nchar(thisLine))
			} else {
				thisStars <- sub( "^ +", "", thisLine)
			}
			starText <- paste( starText, thisStars, sep="")
			curLine <- curLine + 1
			next
		}
		if ( thisLine == "" || substr( thisLine, 1, 4) == "    " || curLine %in% gapLines) {
			curLine <- curLine + 1
			gptr <- 1
			next
		}

		thisGene <- sub( "(^.+)( +)(.+$)","\\1", thisLine)
		thisSeq <- sub( "(^.+)( +)(.+$)","\\3", thisLine)
		thisGene <- gsub( " ", "", thisGene)
		thisSeq <- gsub( " ", "", thisSeq)
		if (gptr > NG) { 
			NG <- gptr
			big[NG] <- gnames[NG] <- ""
		}
		gnames[gptr] <- thisGene
		big[gptr] <- paste( big[gptr], thisSeq, sep="")
		gptr <- gptr + 1
		curLine <- curLine + 1
		starHeadPtr <- c( starHeadPtr, regexpr( thisSeq, thisLine, fixed=T))
	}
	allLens <- nchar(big)
	if (verbose) cat( "\n", NG, gnames, "\nLen: ", allLens)

	NB <- max( allLens)
	bigM <- matrix("", nrow=NG, ncol=NB)
	for ( i in 1:NG) {
		thisLen <- allLens[i]
		bigM[ i, 1:thisLen] <- strsplit( big[i], split="")[[1]][1:thisLen]
	}
	rownames(bigM) <- gnames
	colnames(bigM) <- 1:NB

	starText <- strsplit( starText, split="")[[1]]
	tbl <- table( starText)
	whichStar <- match( "*", names(tbl), nomatch=0)
	pctConserved <- 0
	if (whichStar > 0) pctConserved <- tbl[whichStar] / sum(tbl)
	if (verbose) cat( "\nPct_Conserved: ", as.percent( pctConserved))

	return( list( "alignment"=bigM, "consensus"=starText, "pct.conserved"=pctConserved))
}
