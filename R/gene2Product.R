# gene2Product.R

`gene2Product` <- function( gNames) {

	# default behavior is to return empty character strings
	out <- rep( "", times=length( gNames))

	# make sure set up
	geneMap <- getCurrentGeneMap()
	if ( nrow( geneMap) < 1) {
		cat( "\nWarning:  no gene map is current...")
		return( out)
	}

	# find those genes, and extract out the Product
	where <- base::match( gNames, geneMap$GENE_ID, nomatch=0)
	found <- where > 0
	out[ found] <-  geneMap$PRODUCT[ where]

	# if not yet found, try by common name too
	try2 <- which( !found)
	where <- base::match( gNames[try2], geneMap$NAME, nomatch=0)
	out[ try2[ where > 0]] <- geneMap$PRODUCT[ where]

	if ( "ORIG_ID" %in% colnames(geneMap)) {
		found <- (out != "")
		try3 <- which( !found)
		where <- base::match( gNames[try3], geneMap$ORIG_ID, nomatch=0)
		out[ try3[ where > 0]] <- geneMap$PRODUCT[ where]
	}

	return(out)
}


`gene2ProductAllSpecies` <- function( gNames, hints=NULL) {

	# default behavior is to return empty character strings
	out <- rep( "", times=length( gNames))
	NG <- length(gNames)
	gFac <- factor( gNames)
	gLevels <- levels(gFac)
	NL <- nlevels(gFac)
	gLprod <- rep( "", NL)
	gPtr <- tapply( 1:NG, gFac, FUN=NULL)

	# do for all species
	saveSpecies <- getCurrentSpecies()
	allSpecies <- getCurrentTargetSpecies()

	for ( spec in allSpecies) {
		setCurrentSpecies( spec)
		geneMap <- getCurrentGeneMap()

		# find those genes, and extract out the Product

		# pass 1: perfect hits...
		where <- base::match( gLevels, geneMap$GENE_ID, nomatch=0)
		found <- (where > 0)
		gLprod[ found] <-  geneMap$PRODUCT[ where]
		if ( all( gLprod != "")) break

		if ( "ORIG_ID" %in% colnames(geneMap)) {
			try2 <- which( !found)
			where <- base::match( gLevels[try2], geneMap$ORIG_ID, nomatch=0)
			gLprod[ try2[ where > 0]] <-  geneMap$PRODUCT[ where]
		}

		# pass 2:  human only short names
		if ( spec == "Hs_grc") {
			gNamesTry <- shortGeneName( gLevels,keep=1)
			gNamesMap <- shortGeneName( geneMap$GENE_ID,keep=1)
			where <- base::match( gNamesTry, gNamesMap, nomatch=0)
			found2 <- (where > 0)
			gLprod[ found2] <-  geneMap$PRODUCT[ where]
			if ( all( gLprod != "")) break
			found <- ( found | found2)
		}

		# pass 3: partial matching
		try3 <- which( !found)
		if ( length( try3) > 0) {
			where <- pmatch( gLevels[try3], geneMap$GENE_ID, nomatch=0, duplicates.ok=TRUE)
			found3 <- (where > 0)
			gLprod[ try3[ found3]] <-  geneMap$PRODUCT[ where]
			if ( all( gLprod != "")) break
			found[ try3[ found3]] <- TRUE
		}

		# if not yet found, try by common name too
		try4 <- which( !found)
		if ( length( try4) > 0) {
			where <- base::match( gLevels[try4], geneMap$NAME, nomatch=0)
			gLprod[ try4[ where > 0]] <- geneMap$PRODUCT[ where]
			if ( all( gLprod != "")) break
		}
	}

	# those still blank could be VSA's?
	stillBlank <- which( gLprod == "")
	if ( length( stillBlank) > 0) {
		outvsa <- vsaGeneProduct( gLevels[stillBlank])
		gLprod[ stillBlank] <- outvsa
		stillBlank <- which( gLprod == "")
	}

	# push from the factor back to the real things
	out <- gLprod[ gPtr]

	# if given a set of hints, use that to fill in 'still blank' ones
	if ( length( stillBlank) > 0) {
	    if ( !is.null(hints) && (length(hints) == length( gNames))) {
		stillBlank <- which( out == "")
		if ( length( stillBlank) > 0) {
			out[ stillBlank] <- hints[ stillBlank]
		}
	    }
	}

	setCurrentSpecies( saveSpecies)

	return(out)
}

