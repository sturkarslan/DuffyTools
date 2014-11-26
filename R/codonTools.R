# codonTools.R


`codon.defaults` <- function() {
	# the data object is called 'codonMap'
	data( CodonMap, envir=environment())
	CodonEnv[[ "CodonMap"]] <- codonMap
	if ( nrow( codonMap) != 64) warning("CodonMap:  bad amino acid lookup table")
}


`getCodonMap` <- function() return( CodonEnv[[ "CodonMap"]])


# turn DNA to AA fast
DNAtoAA.fast <- function( dna, clipAtStop=TRUE, readingFrames=1:6) {

	dnaV <- strsplit( dna, split="")[[1]]
	nc <- length( dnaV)
	readingFrames <- intersect( 1:6, readingFrames)
	if ( any( readingFrames %in% 4:6)) {
		dnaRV <- myReverseComplement( dna, as.vector=TRUE)
	}
	nFrames <- length(readingFrames)
	out <- rep( "", times=nFrames)
	names( out) <- readingFrames
	if ( nc < 3) return( out)

	#myPaste <- base::paste

	for ( i in 1:length(readingFrames)) {
		k <- readingFrames[i]
		if ( k <= 3) {
			aaV <- quickDNAtoAA( dnaV[ k:nc])
		} else {
			aaV <- quickDNAtoAA( dnaRV[ (k-3):nc])
		}

		if ( clipAtStop) {
			where <- match( STOP_CODON, aaV, nomatch=0)
			if ( where > 0) length(aaV) <- where
		}
		out[i] <- paste( aaV, collapse="")
	}

	return( out)
}


#  turn a DNA sequence into protein ( all three reading frames)
DNAtoAA <- function( dna, strand="+", clipAtStop=FALSE) {

	dna <- base::toupper( dna)
	if ( base::nchar(dna) < 3) return("")

	out <- vector( mode="character", length=3)
	desc <- "F"

	if ( strand == "-") {
		newdna <- myReverseComplement( dna)
		dna <- newdna
		desc <- "R"
	}

	# set up to fill all three reading frames as we go along
	codonMap <- getCodonMap()

	frame <- 1
	last <- base::nchar(dna) - 2
	for ( i in 1:last) {
		triplet <- base::substr( dna, i, (i+2))
		who <- base::which( codonMap$DNA == triplet)
		if ( length(who) > 0) {
			aa <- codonMap$AA[ who[1]]
		} else {
			aa <- "?"
		}
		# append to that frame
		out[ frame] <- base::paste( out[ frame], aa, sep="")
		frame <- frame + 1
		if ( frame > 3) frame <- 1
	}

	# do we truncate peptides based on stop codon
	if (clipAtStop) {
		stopSet <- regexpr( STOP_CODON, out, fixed=TRUE)
		out <- ifelse( stopSet > 0, base::substr( out, 1, (stopSet-1)), out)
	}

	names( out) <- base::paste( desc, 1:3, sep="")
	return( out)
}


quickDNAtoAA <- function( dnaVec) {

	# faster simple case:  'dna' is a vector, in frame, forward strand.
	Ndna <- length(dnaVec)
	if ( Ndna < 3) return("")

	dna <- base::toupper( dnaVec)
	codonMap <- getCodonMap()

	beg <- seq.int( 1, Ndna, 3)
	Naa <- length(beg)
	end <- beg + 2
	if (end[Naa] > Ndna) {
		Naa <- Naa - 1
		length(beg) <- length(end) <- Naa
	}
	out <- rep.int( "?", Naa)

	#myPaste <- base::paste

	triplets <- mapply( beg, end, FUN=function(b,e) { return( paste( dna[b:e], collapse=""))},
			USE.NAMES=FALSE)
	where <- base::match( triplets, codonMap$DNA, nomatch=0)
	out[ where > 0] <- codonMap$AA[ where]
	return( out)
}


myReverseComplement <- function( dna, as.vector=FALSE) {

	if ( length(dna) > 1) {
		warning( "'reverseComplement' requires a single charater string...dropping some.")
		dna <- dna[1]
	}
	nc <- base::nchar(dna)
	dnaV <- base::unlist( strsplit( dna, ""))
	newV <- rev( dnaV)
	back <- nc:1
	newV[ dnaV[ back] == "A"] <- "T"
	newV[ dnaV[ back] == "C"] <- "G"
	newV[ dnaV[ back] == "G"] <- "C"
	newV[ dnaV[ back] == "T"] <- "A"
	newV[ dnaV[ back] == "N"] <- "N"
	if (as.vector) return( newV)
	return( base::paste( newV, collapse=""))
}


myReverse <- function( dna, as.vector=FALSE) {

	if ( length(dna) > 1) {
		warning( "'reverse' requires a single charater string...dropping some.")
		dna <- dna[1]
	}
	dnaV <- base::unlist( strsplit( dna, ""))
	newV <- rev( dnaV)
	if (as.vector) return( newV)
	return( base::paste( newV, collapse=""))
}


fastAAreadingFrame <- function( peptideTriple, protein, AAprotein, max.mismatch=3) {

	# can we find the best reading frame, given a 'target' protein?
	require( Biostrings)

	# try a fast way first, see if one exact match
	hits <- vector( length=3)
	hits[1] <- regexpr( peptideTriple[1], protein, fixed=TRUE)
	hits[2] <- regexpr( peptideTriple[1], protein, fixed=TRUE)
	hits[3] <- regexpr( peptideTriple[1], protein, fixed=TRUE)
	if ( sum( hits > 0) == 1) {
		best <- which.max( hits)
		bestPep <- peptideTriple[best]
		bestFirstAA <- hits[best]
		return( list( "Position"=bestFirstAA, "Length"=base::nchar(bestPep),
			"Mismatches"=0))
	}

	# see how long the AA chain is...
	firstStopCodon <- regexpr( STOP_CODON, peptideTriple, fixed=TRUE)
	lenAA <- ifelse( firstStopCodon > 0, firstStopCodon - 1, base::nchar( peptideTriple))

	# count up the size of the matches, with various mismatch threshold...
	nMismatch <- 0
	repeat {
		scores <- firstAAinProtein <- lenMatch <- nMisMatch <- rep( 0, times=3)
		for (j in 1:3) {
			pep <- peptideTriple[j]
			len <- lenAA[j]
			if ( len > 0) {
				mp <- matchPattern( pep, AAprotein, max.mismatch=nMismatch)
				if (length(mp) > 0) {
					lens <- nchar( mp)
					bestView <- which.max( lens)
					bestLen <- lens[bestView]
					scores[j] <- bestLen / width(mp)[bestView]
					firstAAinProtein[j] <- start(mp)[bestView]
					lenMatch[j] <- bestLen
					nMisMatch[j] <- nMismatch
				}
			}
		}
		# perfect result is one peptide exactly 1.0, and two at 0.0
		who1 <- (scores == 1.0)
		if ( sum(who1) == 1) {
			curBest <- which.max( scores)
			break
		}
		# not quite perfect,
		curBest <- which.max( scores)
		nMismatch <- nMismatch + 1
		if (nMismatch > max.mismatch) break
	}
	bestPep <- peptideTriple[curBest]
	bestFirstAA <- firstAAinProtein[ curBest]
	return( list( "Position"=bestFirstAA, "Length"=lenMatch[ curBest], 
			"Mismatches"=nMisMatch[curBest]))
}


bestAAreadingFrame <- function( peptideTriple, protein, max.mismatch=3) {

	# can we find the best reading frame, given a 'target' protein?
	require( Biostrings)

	# try a fast way first, see if one exact match
	hits <- vector( length=3)
	hits[1] <- regexpr( peptideTriple[1], protein, fixed=TRUE)
	hits[2] <- regexpr( peptideTriple[1], protein, fixed=TRUE)
	hits[3] <- regexpr( peptideTriple[1], protein, fixed=TRUE)
	if ( sum( hits > 0) == 1) {
		best <- which.max( hits)
		bestPep <- peptideTriple[best]
		bestFirstAA <- hits[best]
		return( list( "Peptide"=bestPep, "Position"=bestFirstAA, "Length"=base::nchar(bestPep),
			"Mismatches"=0))
	}

	# see how long the AA chain is...
	firstStopCodon <- regexpr( STOP_CODON, peptideTriple, fixed=TRUE)
	lenAA <- ifelse( firstStopCodon > 0, firstStopCodon - 1, base::nchar( peptideTriple))

	# count up the size of the matches, with various mismatch threshold...
	nMismatch <- 0
	repeat {
		scores <- rep( 0, times=3)
		firstAAinProtein <- rep( 0, times=3)
		lenMatch <- rep( 0, times=3)
		nMisMatch <- rep( 0, times=3)
		for (j in 1:3) {
			pep <- peptideTriple[j]
			len <- lenAA[j]
			if ( len > 0) {
				mp <<- matchPattern( pep, protein, max.mismatch=nMismatch)
				if (length(mp) > 0) {
					bestView <- which.max( nchar( mp))
					scores[j] <- nchar(mp)[bestView] / width(mp)[bestView]
					firstAAinProtein[j] <- start(mp)[bestView]
					lenMatch[j] <- nchar(mp)[bestView]
					nMisMatch[j] <- nMismatch
				}
			}
		}
		# perfect result is one peptide exactly 1.0, and two at 0.0
		who1 <- (scores == 1.0)
		if ( sum(who1) == 1) {
			curBest <- which.max( scores)
			break
		}
		# not quite perfect,
		curBest <- which.max( scores)
		nMismatch <- nMismatch + 1
		if (nMismatch > max.mismatch) break
	}
	bestPep <- peptideTriple[curBest]
	bestFirstAA <- firstAAinProtein[ curBest]
	#cat( "\nScores: ", scores, "\tBest ReadingFrame: \t", bestPep, "\tStart location:\t", bestFirstAA)
	return( list( "Peptide"=bestPep, "Position"=bestFirstAA, "Length"=lenMatch[ curBest], 
			"Mismatches"=nMisMatch[curBest]))
}


# try to convert DNA to AA given the current annotation info
`convertGenomicBasesToCodingAminoAcids` <- function( seqID, position, end, strand="+", dnaQuery, genomeDNA,
			geneMap=NULL, exonMap=NULL) {

	# given a vector of DNA bases and the its genomic location  < dnaQuery, position, end >,
	# convert that to the amino acid sequence for the coding strand that covers...
	nBases <- length( dnaQuery)
	outG <- outQ <- rep( "", times=nBases)
	out <- list( "genomic"=outG, "query"=outQ)

	if ( nBases != (end - position + 1)) {
		warning( "bad DNA query string and/or size")
		return(out)
	}

	# get the gene(s) covered
	if ( is.null( geneMap)) geneMap <- getCurrentGeneMap()
	gmap <- subset.data.frame( geneMap, SEQ_ID == seqID & POSITION < end & END > position & REAL_G == TRUE)
	if ( nrow(gmap) < 1) return( out)

	if ( is.null( exonMap)) exonMap <- getCurrentExonMap()
	exonMap <- subset.data.frame( exonMap, GENE_ID %in% gmap$GENE_ID)

	for( ig in 1:nrow( gmap)) {
		thisG <- gmap$GENE_ID[ ig]
		thisStrand <- gmap$STRAND[ ig]
		emap <- subset.data.frame( exonMap, GENE_ID == thisG)
		if ( nrow(emap) < 1) next

		# build the string of coding nucleotides from the forward strand
		forwardDNA <- vector()
		for ( ie in 1:nrow(emap)) {
			thisExon <- genomeDNA[ emap$POSITION[ie] : emap$END[ie] ]
			names( thisExon) <- emap$POSITION[ie] : emap$END[ie] 
			forwardDNA <- base::append( forwardDNA, thisExon)
		}
		codingDNA <- forwardDNA
		if ( thisStrand == "-") {
			tmp <- base::paste( codingDNA, collapse="")
			tmp <- myReverseComplement(tmp)
			tmp <- strsplit( tmp, split="")[[1]]
			codingDNA <- tmp
			names( codingDNA) <- rev( names( forwardDNA))
		}

		# force a trim to multiple of 3 bases
		bigN <- floor( length(codingDNA)/3) * 3
		nCodingBases <- length( codingDNA) <- bigN
		if ( nCodingBases < 3) next

		# convert to AA, ( we know we are in Frame, so just keep the first reading frame)
		thisAA <- quickDNAtoAA( codingDNA)

		# build this back into a vector of single letters, with the AA at the center of its 3 bases
		newstr <- rep( "", times=nCodingBases)
		newstr[ seq( 2, nCodingBases, by=3)] <- thisAA
		names( newstr) <- names( codingDNA)

		# lastly, we can put these AA back into the correct location of the output genomic string of text
		for ( ic in 1:nCodingBases) {
			if ( newstr[ic] == "") next
			genLocation <- as.integer( names( newstr)[ic])
			outLocation <- genLocation - position + 1
			if ( outLocation < 1 || outLocation > nBases) next
			outG[ outLocation] <- newstr[ic]
		}

		# now put the query bases in, and do it again
		outQ <- outG
		codingDNA <- forwardDNA
		codingLocs <- as.integer( names( codingDNA))
		queryLocs <- position : end
		where <- base::match( codingLocs, queryLocs, nomatch=0)
		codingHits <- which( where > 0)
		if ( length( codingHits) > 0) {
			newstr2 <- newstr
			codingDNA[ codingHits] <- dnaQuery[ where]
			# with the possibility of Indels, these query bases may not be singletons any more
			isDiff <- which( codingDNA != forwardDNA)
			if ( length(isDiff)) {
				baseLen <- nchar( codingDNA)
				# try to estimate how many 'post-indel' AA for each 'pre-indel' codon, pad the ends for edge cases
				baseCumSum <- cumsum( c( 1, baseLen, 1))
				snpDNA <- base::paste( codingDNA, collapse="")
				if ( thisStrand == "-") {
					snpDNA <- myReverseComplement(snpDNA)
					baseCumSum <- cumsum( c( 1, rev(baseLen), 1))
				}
				snpDNA <- strsplit( snpDNA, split="")[[1]]
				if ( thisStrand == "-") {
					names( snpDNA) <- rep( rev( names( forwardDNA)), times=baseLen)
				} else {
					names( snpDNA) <- rep( names( forwardDNA), times=baseLen)
				}
				snpAA <- quickDNAtoAA( snpDNA)
				# this length may be different due to indels, step along by hand...
				newNbases <- min( nCodingBases, length(snpDNA))
				newNaa <- floor( newNbases/3)
				nAAused <- 0
				for ( j in seq( 2, newNbases, by=3)) {
					# account for the padding around the cum sum vector
					ndna <- diff( baseCumSum[ c(j-1,j+2)])
					naa <- round( ndna/3)
					if (naa == 1) {
						newstr2[ j] <- snpAA[nAAused+1]
						nAAused <- nAAused + 1
					} else if ( naa < 1) {
						newstr2[ j] <- ""
					} else {
						newstr2[ j] <- paste( snpAA[ (nAAused+1):(nAAused+naa)], collapse="")
						if ( thisStrand == "-") newstr2[j] <- myReverse( newstr2[j])
						nAAused <- nAAused + naa
					}
				}
				whodiff <- which( newstr2 != newstr)
				for ( ic in whodiff) {
					if ( newstr2[ic] == "") next
					genLocation <- as.integer( names( newstr2)[ic])
					outLocation <- genLocation - position + 1
					if ( outLocation < 1 || outLocation > nBases) next
					outQ[ outLocation] <- newstr2[ic]
				}
			} else {
				# no Indels, so do it fast and like before...
				if ( thisStrand == "-") {
					tmp <- base::paste( codingDNA, collapse="")
					tmp <- myReverseComplement(tmp)
					tmp <- strsplit( tmp, split="")[[1]]
					codingDNA <- tmp
					names( codingDNA) <- rev( names( forwardDNA))
				}
				thisAA <- quickDNAtoAA( codingDNA)
				newstr2[ seq( 2, nCodingBases, by=3)] <- thisAA
				whodiff <- which( newstr2 != newstr)
				for ( ic in whodiff) {
					if ( newstr2[ic] == "") next
					genLocation <- as.integer( names( newstr2)[ic])
					outLocation <- genLocation - position + 1
					if ( outLocation < 1 || outLocation > nBases) next
					outQ[ outLocation] <- newstr2[ic]
				}
			}
		}
	}

	out <- list( "genomic"=outG, "query"=outQ)
	return( out)
}


# try to convert AA to DNA given the current annotation info
`convertAApositionToGenomicDNAposition` <- function( geneID, AAposition, AAlength) {

	outPos <- outEnd <- 0
	out <- list( "SEQ_POSITION"=outPos, "SEQ_END"=outEnd)

	emap <- subset.data.frame( getCurrentExonMap(), GENE_ID == geneID)
	if ( nrow( emap) < 1) return( out)

	# make a little band of relative DNA positions
	vNow <- vector()
	for ( j in 1:nrow(emap)) {
		thisBeg <- emap$POSITION[j]
		thisEnd <- emap$END[j]
		sml <- thisBeg:thisEnd
		vNow <- base::append( vNow, sml)
	}
	names( vNow) <- 1:length(vNow)
	if ( emap$STRAND[1] == "-") names( vNow) <- rev( 1:length(vNow))

	firstAAbase <- (AAposition-1) * 3 + 1
	lastAAbase <- (AAposition+AAlength-1) * 3 

	where <- base::match( firstAAbase, names(vNow), nomatch=0)
	if ( where > 0) outPos <- vNow[ where]
	where <- base::match( lastAAbase, names(vNow), nomatch=0)
	if ( where > 0) outEnd <- vNow[ where]
	if ( outPos > outEnd) { tmp <- outPos; outPos <- outEnd; outEnd <- tmp }
	names(outPos) <- names(outEnd) <- ""

	out <- list( "SEQ_POSITION"=outPos, "SEQ_END"=outEnd)
	return( out)
}


# try to convert genomic DNA to coding DNA given the current annotation info
`convertGenomicDNAtoCodingDNA` <- function( geneID, genomicDNA=NULL) {

	outDNA <- ""

	emap <- subset.data.frame( getCurrentExonMap(), GENE_ID == geneID)
	if ( nrow( emap) < 1) return( out)
	if ( is.null( genomicDNA)) {
		cat( "\nRequired 'genomicDNA' argument is missing...")
		return( out)
	}

	seqID <- emap$SEQ_ID[1]

	# make a little band of absolute DNA positions
	vNow <- vector()
	for ( j in 1:nrow(emap)) {
		thisBeg <- emap$POSITION[j]
		thisEnd <- emap$END[j]
		sml <- thisBeg:thisEnd
		vNow <- base::append( vNow, sml)
	}

	# get that chunk of genomic DNA
	beg <- min( vNow)
	end <- max( vNow)
	myDNA <- strsplit( substr( genomicDNA, beg, end), split="")[[1]]

	# convert to relative positions
	vNow <- vNow - beg + 1
	myDNA <- myDNA[ vNow]
	outDNA <- base::paste( myDNA, collapse="")

	# flip if reverse strand
	if ( emap$STRAND[1] == "-") outDNA <- myReverseComplement( outDNA)

	return( outDNA)
}


`buildCodonFreqMap` <- function( AAfasta, speciesID="Hs_grc", fraction=1.0, genomicFastaFile=NULL) {

	codonTable <- measureCodonFreq( AAfasta=AAfasta, speciesID=speciesID, fraction=fraction, 
					genomicFastaFile=genomicFastaFile)
	codonMap <- calculateCodonScores( codonTable)

	return( list( "table"=codonTable, "map"=codonMap))
}


`measureCodonFreq` <- function( AAfasta, speciesID="Hs_grc", fraction=1.0, 
				GImap=NULL, genomicFastaFile=NULL) {

	setCurrentSpecies( speciesID)
	geneMap <- getCurrentGeneMap()
	geneMap <- subset.data.frame( geneMap, REAL_G == TRUE)

	# given a pair of matching fasta files, one amino acids and one cDNA,
	# build the table of most frequent codons for each AA

	aaList <- loadFasta( AAfasta)
	aaNames <- aaList$desc
	aaNames <- sub( "(gi.+ref\\|)(NP_.+)(\\.[1-9]\\|.+$)", "\\2", aaNames)
	aaSeqs <- aaList$seq
	cat( "\n  N_Proteins:  ", length( aaNames), head( aaNames))

	# convert to GeneIDs now...
	if ( is.null( GImap)) {
		aaGeneIDs <- aaNames
	} else {
		where <- base::match( aaNames, GImap$protAcc, nomatch=0)
		aaGeneIDs <- rep( "", times=length(where))
		thisGeneName <- GImap$name[where]
		thisGInumber <- GImap$locusLinkId[where]
		allGeneID <- base::paste( thisGeneName, thisGInumber, sep=":GI")
		aaGeneIDs[ where > 0] <- allGeneID
	}
	aaGeneIDs <- intersect( aaGeneIDs, geneMap$GENE_ID)
	cat( "\n  N_Matching_Genes:   ", length( aaGeneIDs), head( aaGeneIDs))

	cat( "\nLooking up genes in geneMap...\n")
	aaGenePtrs <- sapply( aaGeneIDs, function(x) {
				if ( x == "") return( 0)
				where <- grep( x, geneMap$GENE_ID, fixed=T)[1]
				if ( is.na(where)) return( 0)
				return( where)
			})
	cat( "\nN_found:  ", sum( aaGenePtrs > 0), "  N_missing:  ", sum( aaGenePtrs == 0))


	codonMap <- getCodonMap()
	codonTable <- vector()

	visitOrder <- base::order( aaGenePtrs)

	cat( "\n\nVisiting all named proteins")
	curSeqID <- ""
	genomicDNA <- ""

	for (i in visitOrder) {

		if ( i == 0) next
		thisname <- aaNames[i]
		thisGeneID <- aaGeneIDs[i]
		thisAA <- aaSeqs[i]
		thisGenePtr <- aaGenePtrs[i]
		if ( thisGenePtr < 1) next
		geneID <- geneMap$GENE_ID[ thisGenePtr]
		seqID <- geneMap$SEQ_ID[ thisGenePtr]
		if ( seqID != curSeqID) {
			genomicDNA <- getFastaSeqFromFilePath( filePath=genomicFastaFile, seqID=seqID)
			curSeqID <- seqID
		}
		thisDNA <- convertGenomicDNAtoCodingDNA( geneID, genomicDNA)

		# for this AA and DNA, see what codon is for each AA
		allAA <- strsplit( thisAA, split="")[[1]]
		allDNA <- strsplit( thisDNA, split="")[[1]]

		allCODONs <- sapply( seq( 1, length(allDNA), by=3), function(x) {
					base::paste( allDNA[ x:(x+2)], collapse="")
				})
		nUse <- round( fraction * min( length( allAA), length(allCODONs)))

		smallTable <- base::table( base::paste( allAA[1:nUse], allCODONs[1:nUse], sep=":"))
		if ( length( codonTable) > 0) {
			codonTable <- mergeTables( codonTable, smallTable)
		} else {
			codonTable <- smallTable
		}

		if ( i %% 200 == 0) {
			cat( "\n", i, thisname, "\n")
			print( head( ibase::sort( codonTable, decreasing=T)))
		}
	}
	cat( "\n")

	return( codonTable)
}


`calculateCodonScores` <- function( codonTable) {

	hasN <- grep( ":.*N", names( codonTable))
	if ( length(hasN) > 0) codonTable <- codonTable[ -hasN]

	# now get the top codons for each AA
	AAout <- sub( ":.+", "", names( codonTable))

	allAA <- allCodons <- allCounts <- allPercents <- allPct1 <- allPct2 <- allPct3 <- vector()
	nCodons <- 0
	
	tapply( 1:length(AAout), INDEX=factor( AAout), FUN=function(x) {

			# given the set of rows in 'codonTable' that all share one AA
			# turn the relative abundance of each base into a phred score
			m <- matrix( 0, nrow=4, ncol=3)
			rownames(m) <- c( "A","C","G","T")
			sumCnts <- sum( codonTable[ x])
			myCodons <- myX <- myCnts <- vector()
			myBases <- vector( mode="list")
			nNow <- 0
			for ( i in x) {
				thisCnt <- codonTable[i]
				# keep all that are more than 5% of the time
				if ( thisCnt < sumCnts * 0.05) next
				thisCodon <- sub( ".:", "", names( codonTable)[i])
				theseBases <- strsplit( thisCodon, split="")[[1]]
				if ( ! all( theseBases %in% rownames(m))) next
				m[ theseBases[1], 1] <-  m[ theseBases[1], 1] + thisCnt
				m[ theseBases[2], 2] <-  m[ theseBases[2], 2] + thisCnt
				m[ theseBases[3], 3] <-  m[ theseBases[3], 3] + thisCnt
				nNow <- nNow + 1
				myX[nNow] <- i
				myCodons[nNow] <- thisCodon
				myBases[[nNow]] <- theseBases
				myCnts[nNow] <- thisCnt
			}
			totalCnts <- apply( m, MARGIN=2, FUN=sum)
			visitOrd <- base::order( myCnts, decreasing=T)
			for ( j in 1:length(myX)) {
				i <- visitOrd[j]
				thisCodon <- myCodons[i]
				thisCnt <- myCnts[i]
				thisBases <- myBases[[i]]
				scores <- c( 0, 0, 0)
				for ( j in 1:3) {
					if ( ! thisBases[j] %in% rownames(m)) next
					scores[j] <- round( m[ thisBases[j], j] * 100 / totalCnts[j])
				}
				nCodons <<- nCodons + 1
				allAA[nCodons] <<- AAout[myX[1]]
				allCodons[nCodons] <<- thisCodon
				allCounts[nCodons] <<- thisCnt
				allPercents[nCodons] <<- thisCnt * 100 / sumCnts
				allPct1[nCodons] <<- scores[1]
				allPct2[nCodons] <<- scores[2]
				allPct3[nCodons] <<- scores[3]
			}
			return()
		})

	# turn the list into a Nx4 matrix
	basePcts <- data.frame( "Pct1"=allPct1, "Pct2"=allPct2, "Pct3"=allPct3, stringsAsFactors=F)
	rownames(basePcts) <- 1:nrow(basePcts)
	phredIntScore <- apply( basePcts, MARGIN=1, function(x) {
					phred <- round( x * 40 / 100)
					return( base::paste( phred, collapse=" "))
				})

	phredAsciiScore <- sapply( phredIntScore, solexaToPhred, scoreType="Phred33")

	out <- data.frame( "AA"=allAA, "Codon"=allCodons, "Count"=allCounts, "Percent"=allPercents,
				basePcts, "PhredIntegers"=phredIntScore, "PhredString"=phredAsciiScore)
	return( out)
}


`test.codonTools` <- function() {
	checkEquals( strsplit(DNAtoAA.fast( "ATGTAG"),split="")[[1]], 
			quickDNAtoAA( strsplit("ATGTAG", split="")[[1]]))
	checkEquals( DNAtoAA( "ATGTAG"), c("F1"="M*", "F2"="C", "F3"="V"))
}
