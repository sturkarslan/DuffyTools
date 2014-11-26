# orthologTools.R


`ortholog` <- function( genes, from="PF3D7", to="PY17X") {

	noOrtho <- ""
	gOut <- rep( noOrtho, times=length(genes))

	# do we need to reload ?
	if ( ! exists( "OrthoTable", envir=OrthoEnv)) {
		orthoTable <- NULL
		toLoad <- "Pf.OrthologTable"
		data( list=list( toLoad), envir=environment())
		if ( is.null( orthoTable)) {
			cat( "\nFailed to load Ortholog Table:  ", toLoad)
			OrthoEnv[[ "OrthoTable"]] <- NULL
			return( gOut)
		}
		
		# prep it a bit...
		OrthoEnv[[ "OrthoTable"]] <- orthoTable
		OrthoEnv[[ "OrthoStrains"]] <- colnames(orthoTable)
	}

	orthoTable <- OrthoEnv[[ "OrthoTable"]]
	orthoStrains <- OrthoEnv[[ "OrthoStrains"]]
	ifrom <- match( toupper(from), orthoStrains, nomatch=0)
	if ( ifrom == 0) {
		cat( "\n'from' strain not in Ortholog table: ", from)
		cat( "\nKnown strains: ", orthoStrains)
		return( gOut)
	}
	ito <- match( toupper(to), orthoStrains, nomatch=0)
	if ( ito == 0) {
		cat( "\n'to' strain not in Ortholog table: ", to)
		cat( "\nKnown strains: ", orthoStrains)
		return( gOut)
	}

	where <- match( gsub( " ", "", genes), orthoTable[ , ifrom], nomatch=0)
	gOut[ where > 0] <- orthoTable[ where, ito]

	return( gOut)
}

