# fastaTools.R  --  collection of FASTA file manipulation routines


`loadFasta` <- function( file="file.fasta", verbose=TRUE, short.desc=TRUE) {

	require( Biostrings)

	file <- allowCompressedFileName( file)
	if (verbose) cat( "\nLoading Fasta file: ", file, "...")

	fa <- readBStringSet( file)
	seqs <- as.character(fa, use.names=TRUE)
	nams <- names(seqs)
	names(seqs) <- NULL

	# for consistency with other tools, trim the descriptor after the first blank
	if ( short.desc) {
		nams <- sub( " .+", "", nams)
	}

	return( list( "desc"=nams, "seq"=seqs))
}


`as.Fasta.data.frame` <- function( fasta) {

	return( data.frame( "desc"=fasta$desc, "seq"=fasta$seq, stringsAsFactors=FALSE))
}


`as.Fasta` <- function( desc, seq) {

	if ( length( desc) != length( seq)) stop( "as.Fasta:  unequal length arguments")

	return( list( "desc"=desc, "seq"=seq))
}


`writeFasta` <- function( fasta, file=NULL, line.width=80) {

	if ( is.null( file)) stop( "writeFasta:  required 'file' argument is missing")
	writeLines( as.text.Fasta( fasta, line.width=line.width), con=file, sep="\n")
}


`writeLongFasta` <- function( desc, seq, file=NULL) {

	if ( is.null( file)) stop( "writeLongFasta:  required 'file' argument is missing")
	writeLines( base::paste( ">", desc, "\n", seq, sep=""), con=file, sep="\n")
}


`as.text.Fasta` <- function( fasta, line.width=80) {

	if ( ! is.list( fasta)) return("")

	N <- length(fasta$desc)
	if ( is.null(N) || N < 1) return("")

	out <- sapply( 1:N, function(i) {
		nam <- fasta$desc[i]
		seq <- fasta$seq[i]
		base::paste( ">", fasta$desc[i], "\n", wrap.text( fasta$seq[i], line.width=line.width), sep="")
	})
	return( out)
}


wrap.text <- function( txt, line.width=60) {

	# re-format fasta text to be line wrapped to a fixed width
	N <- length(txt)
	if ( is.null(txt) || N < 1) return("")

	out <- sapply( 1:N, function(i) {
			newtxt <- ""
			oldtxt <- txt[i]
			nch <- base::nchar( oldtxt)
			if ( nch < 1) return("")
			smltxt <- sapply( seq( 1, nch, by=line.width), function(j) {
					last <- min( (j+line.width-1), nch)
					base::substr( oldtxt, j, last)
				})
			paste( smltxt, collapse="\n")
		})
	return( out)
}


# smart fasta file lookup...
FastaFilePathEnv <- new.env( parent=emptyenv())
assign( "currentFastaFile", "", envir=FastaFilePathEnv)
assign( "currentFastaObject", NULL, envir=FastaFilePathEnv)


# get one FASTA sequence, by filename and seqID.  returns a Biostrings object.
`getFastaSeqFromFilePath` <- function( filePath, seqID, verbose=FALSE) {

	# see if we need to read a different file
	alreadyLoaded <- ( filePath == get( "currentFastaFile", envir=FastaFilePathEnv))

	if ( verbose) cat( "\nGetting FASTA seq for: ",seqID)
	if ( ! alreadyLoaded) {
		# we could be given an explicit filename OR a directory
		info <- file.info( filePath)
		if ( any( is.na( info$size))) stop( paste( "getFastaSeqFromFilePath:  file not found:  ", filePath))
		isDirectory <- info$isdir
		if( isDirectory) {
			pathRelative <- TRUE
			files <- dir( filePath)
			# try to find a file that has that seqID as part of its name
			if (verbose) cat( "   trying", length(files), "files in folder.")
			curSpecies <- getCurrentSpecies()
			tryFileName <- paste( seqID, ".fa", sep="")
			hit <- pmatch( tryFileName, files, nomatch=0)
			if ( hit == 0) {
				if (verbose) cat( "   trying prepend of speciesID.")
				tryFileName <- sub( paste( curSpecies,"_",sep=""), "", tryFileName, fixed=TRUE)
				hit <- pmatch( tryFileName, files, nomatch=0)
			}

			if ( hit == 0) {
			    files <- dir( filePath, full.name=T)
			    # last chance:  see if any subfolders have that file
			    myfolders <- files[ file.info( files)$isdir]
			    if ( length( myfolders) > 0) {
				pathRelative <- FALSE
			    	if (verbose) cat( "   trying", length( myfolders), "subfolders.")
			  	morefiles <- vector()
			  	for( f in myfolders) morefiles <- append( morefiles, dir( f, full.name=T))
			  	tryFileName <- paste( seqID, ".fa", sep="")
				#hit <- pmatch( tryFileName, morefiles, nomatch=0)
				hit <- grep( tryFileName, morefiles)
				hit <- if ( length(hit) > 0) hit[1] else 0
				if ( hit == 0) {
				    if (verbose) cat( "   trying prepend of speciesID.")
				    tryFileName <- sub( paste( curSpecies,"_",sep=""), "", tryFileName, fixed=TRUE)
				    #hit <- pmatch( tryFileName, morefiles, nomatch=0)
				    hit <- grep( tryFileName, morefiles)
				    hit <- if ( length(hit) > 0) hit[1] else 0
				}
				files <- morefiles
			    }
			}

			if ( hit == 0) stop( paste( "\nUnable to find FASTA file for:  ", seqID, "  in folder:  ", filePath))

			if (pathRelative) {
				useFile <- file.path( filePath, files[ hit])
			} else {
				useFile <- files[ hit]
			}
		} else {
			useFile <- filePath
		}

		# ok, we have a file, so load it
		assign( "currentFastaObject", loadFasta( useFile), envir=FastaFilePathEnv)
		assign( "currentFastaFile", useFile, envir=FastaFilePathEnv)
	}

	# now look for that seqID
	fasta <- get( "currentFastaObject", envir=FastaFilePathEnv)
	where <- match( seqID, fasta$desc, nomatch=0)
	if ( where > 0) return( fasta$seq[ where])

	# complain and quit if not found
	warning( paste( "Fasta descriptor:  ", seqID, "   not found in Fasta file/path:  ", filePath))
	return( NULL)
}
