# humanIDterms.R - turn human GENE_IDs into two extra columns ENTREZ_ID and common NAME


getHumanIDterms <- function( geneIDs) {

	geneIDs <- as.character( geneIDs)

	# format is {commonname:GInumber:chromosome:location}
	ginum <- sub( "(^.+:GI)([0-9]+)(:?.*$)", "\\2", geneIDs)
	nam <- sub( "(^.+?)(:.*$)", "\\1", geneIDs)

	# verify the Entrez ID is valid...
	suppressWarnings( ginum[ is.na( as.integer( ginum))] <- "" )

	out <- list( "GENE_NAME"=nam, "ENTREZ_ID"=ginum)
	return( out)
}


addHumanIDterms <- function( mydf, idColumn="GENE_ID") {

	if ( ! idColumn %in% colnames(mydf)) {
		cat( "\nHuman GeneID column not found: ", idColumn, "\nFound: ", colnames(mydf))
		return( mydf)
	}

	humanTerms <- getHumanIDterms( mydf[[ idColumn]])

	out <- cbind( as.data.frame( humanTerms), mydf, stringsAsFactors=FALSE)
	return( out)
}
