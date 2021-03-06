# targetTools.R

# utilities to manage the Target Organisms and Indexes inside the DuffyRNAseq package


# quick wrappers to get one 
`getCurrentTarget` <- function() { 
	
	ans <- TargetEnv[[ "CurrentTarget" ]]
	speciesSet <- getCurrentTargetSpecies()
	txt <- sapply( speciesSet, FUN=getSpeciesText)
	txt <- paste( txt, collapse="; ")
	ans$SpeciesText <- txt
	return( ans)
}


`setCurrentTarget` <- function( targetID=NULL, optionsFile=NULL) {

	allTargets <- TargetEnv[[ "AllTargets"]]

	if ( is.null( targetID)) {
		targetID <- getOptionValue( optionsFile, "targetID", notfound="HsPf")
	}
	where <- base::match( targetID, allTargets$TargetID, nomatch=0)
	if ( where == 0) {
		stop( paste( "setCurrentTarget:  Error:  targetID not found.  Given: ", 
				targetID, "Current choices:  ", 
				paste( allTargets$TargetID, collapse=" | ")))
	}

	
	assign( "CurrentTarget", value=allTargets[ where, ],  envir=TargetEnv)

	# force the current species to the first in this target
	speciesSet <- getCurrentTargetSpecies()
	setCurrentSpecies( speciesID=speciesSet[1])

	return( allTargets$TargetID[ where])
}
	


# get the entire Target Set at once
`getAllTargets` <- function() { 

	ans <- TargetEnv[[ "AllTargets"]]
	ans$SpeciesText <- ""
	for ( i in 1:nrow(ans)) {
		speciesSet <- strsplit( ans$SpeciesSet[i], split=",")[[1]]
		txt <- sapply( speciesSet, FUN=getSpeciesText)
		txt <- paste( txt, collapse="; ")
		ans$SpeciesText[i] <- txt
	}
	return(ans)
}


# get the set of speciesIDs for the current target
`getCurrentTargetSpecies` <- function() { 

	thisTarget <- TargetEnv[[ "CurrentTarget" ]]
	if (nrow( thisTarget) < 1) return( vector())
	speciesSet <- strsplit( gsub( " ", "", thisTarget$SpeciesSet), split=",", fixed=T)[[1]]
	return( speciesSet)
}


# get the set of speciesIDs for the current target
`getCurrentTargetFilePrefix` <- function() { 

	thisTarget <- TargetEnv[[ "CurrentTarget" ]]
	if (nrow( thisTarget) < 1) return( vector())
	prefixSet <- strsplit( gsub( " ", "", thisTarget$PrefixSet), split=",", fixed=T)[[1]]
	return( prefixSet)
}


# add an explicit target to the targets envirn
`addTarget` <- function( targetID, speciesSet=targetID, prefixSet=speciesSet, mapset.path=NULL) {

	speciesVec <- strsplit( gsub( " ", "", speciesSet), split=",", fixed=T)[[1]]
	prefixVec <- strsplit( gsub( " ", "", prefixSet), split=",", fixed=T)[[1]]
	if ( length(speciesVec) != length( prefixVec)) {
		cat( "\nNeed the same number of species IDs and file prefixes")
		cat( "\nSpecies:  ", speciesVec)
		cat( "\nPrefixes: ", prefixVec)
		stop()
	}

	# make a one line data.frame for this given target
	thisTarget <- data.frame( targetID, speciesSet, prefixSet, stringsAsFactors=FALSE)
	colnames( thisTarget) <- TARGET_NAMES

	# get the currently known targets and see if its already in
	allTargets <- TargetEnv[[ "AllTargets"]]
	who <- (-1)
	if (nrow( allTargets) > 0) {
		who <- base::match( targetID, allTargets$TargetID, nomatch=0)
	}

	# either overwrite, append, or start fresh
	if ( who > 0) {
		allTargets[ who, ] <- thisTarget[ 1, ]
		myRow <- who
	} else if ( who == 0) {
		allTargets <- rbind( allTargets, thisTarget)
		myRow <- nrow( allTargets)
	} else {
		allTargets <- thisTarget
		myRow <- 1
	}
	rownames( allTargets) <- 1:nrow(allTargets)

	# store the complete set, and the current one
	assign( "AllTargets", value=allTargets, envir=TargetEnv)
	assign( "CurrentTarget", value=allTargets[myRow, ], envir=TargetEnv)

	# force the load of these mapsets
	for ( i in 1:length(speciesVec)) {
		spec <- speciesVec[i]
		if ( spec %in% getAllSpecies()) next
		cat( "  ", spec, "..", sep="")
		mapsetName <- paste( prefixVec[i], "MapSet", sep=".")
		mapSet <- findAndLoadMapSet( mapsetName, mapset.path=mapset.path)
		if ( is.null( mapSet)) {
			cat( "\nFailed to find/load MapSet data...  Species: ", spec, "\tPrefix: ", prefixVec[i], "\n")
		}
	}

	setCurrentSpecies( speciesVec[1])
	return( targetID)
}



# reset to factory defaults
`target.defaults` <- function() {

	# load the package mapSets, and stuff them into the mapSet environment
	assign( "AllTargets", value=data.frame(), envir=TargetEnv)
	assign( "CurrentTarget", value=data.frame(), envir=TargetEnv)

	cat( "\nLoading Target Species: ")
	addTarget( "Pf3D7", "Pf3D7", "Pf")
	addTarget( "Hs_grc", "Hs_grc", "Hs")
	addTarget( "HsPf", "Pf3D7,Hs_grc", "Pf,Hs")

	addTarget( "PCO", "PCO", "Pco")
	addTarget( "PkH", "PkH", "PkH")
	addTarget( "MacMu", "MacMu", "MacMu")
	addTarget( "MacPCO", "PCO,MacMu", "Pco,MacMu")
	addTarget( "MacPkH", "PkH,MacMu", "PkH,MacMu")

	addTarget( "PbANKA", "PbANKA", "Pb")
	addTarget( "PbMmu", "PbANKA,Mmu_grc", "Pb,Mmus")

	addTarget( "Py17X", "Py17X", "Py17X")
	addTarget( "Mmu_grc", "Mmu_grc", "Mmus")
	addTarget( "PyMmu", "Py17X,Mmu_grc", "Py17X,Mmus")
	addTarget( "PyHsMmu", "Py17X,Hs_grc,Mmu_grc", "Py17X,Hs,Mmus")

	addTarget( "PvSal1", "PvSal1", "Pv")
	addTarget( "PvMmu", "PvSal1,Mmu_grc", "Pv,Mmus")
	addTarget( "PvHsMmu", "PvSal1,Hs_grc,Mmu_grc", "Pv,Hs,Mmus")

	addTarget( "Agam", "Agam", "Ag")
	addTarget( "AgPf", "Pf3D7,Agam", "Pf,Ag")

	addTarget( "Ecoli", "Ecoli", "Eco")
	addTarget( "Styphi", "Styphi", "Sty")

	addTarget( "MT_H37", "MT_H37", "MTb")
	addTarget( "HsMTb", "Hs_grc,MT_H37", "Hs,MTb")

	#addTarget( "LmjF", "LmjF", "LmjF")
	#addTarget( "LinJ", "LinJ", "LinJ")
	
	#addTarget( "Msmeg_mc2", "Msmeg_mc2", "Msmeg")
	addTarget( "PfDD2", "PfDD2", "PfDD2")
	addTarget( "PfIT", "PfIT", "PfIT")

	setCurrentTarget( "HsPf")
	cat( "  Done.\n")
	return()
}


`exportTargets` <- function( fileout="DuffyTools.Targets.txt") {

	tmp <- getAllTargets()
	write.table( tmp, file=fileout, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
	cat( "\nWrote targets file:  ", fileout, "\nN_Targets: ", nrow( tmp), "\n")
	return()
}


`importTargets` <- function( filein="DuffyTools.Targets.txt") {


	if ( ! file.exists(filein)) stop( paste( "importTargets:  cannot file targets file",
			"\nTried: ", filein))

	tmp <- read.delim( filein, header=TRUE, as.is=TRUE)

	# re-initialize to empty
	assign( "AllTargets", value=data.frame(), envir=TargetEnv)
	assign( "CurrentTarget", value=data.frame(), envir=TargetEnv)

	for ( i in 1:nrow(tmp)) {
		addTarget( tmp$TargetID[i], tmp$SpeciesSet[i], tmp$PrefixSet[i])
	}

	# use the first as default
	assign( "CurrentTarget", value=tmp$TargetID[1], envir=TargetEnv)
	# force the current species to the first in this target
	speciesSet <- getCurrentTargetSpecies()
	setCurrentSpecies( speciesID=speciesSet[1])

	return( getAllTargets()$TargetID)
}

