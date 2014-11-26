# aliasTools.R


`alias2Gene` <- function( genes) {

	# do we need to reload ?
	curSpecies <- getCurrentSpecies()
	if ( (!exists( "AliasTable", envir=AliasEnv)) || (AliasEnv[[ "AliasSpecies"]] != curSpecies)) {
		AliasTable <- NULL
		toLoad <- paste( getCurrentSpeciesFilePrefix(), "AliasTable", sep=".")
		data( list=list( toLoad), envir=environment())
		if ( is.null( AliasTable)) {
			cat( "\nFailed to load Alias Table:  ", toLoad)
			AliasEnv[[ "AliasTable"]] <- NULL
			AliasEnv[[ "AliasSpecies"]] <- ""
			return( genes)
		}
		
		# prep it a bit...
		AliasTable$Alias <- toupper( AliasTable$Alias)
		AliasTable$Version <- as.numeric( AliasTable$Version)
		ord <- order( AliasTable$Alias, -(AliasTable$Version))
		AliasTable <- AliasTable[ ord, ]
		AliasEnv[[ "AliasTable"]] <- AliasTable
		AliasEnv[[ "AliasSpecies"]] <- curSpecies
	}

	aliasTable <- AliasEnv[[ "AliasTable"]]

	N <- length( genes)
	genesOut <- genes
	versionOut <- rep( 0, times=N)

	repeat {
		where <- base::match( toupper(genesOut), aliasTable$Alias, nomatch=0)
		hasEntry <- which( where > 0)
		if ( length( hasEntry) < 1) break
		newGene <- aliasTable$GeneID[ where]
		newVersion <- aliasTable$Version[ where]
		okToChange <- which( versionOut[hasEntry] < newVersion)
		if ( length( okToChange) < 1) break

		genesOut[ hasEntry[ okToChange]] <- newGene[ okToChange]
		versionOut[ hasEntry[ okToChange]] <- newVersion[ okToChange]
	}

	return( genesOut)
}

