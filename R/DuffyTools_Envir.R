# DuffyTools_Envir.R

# set up and manage the 'global environment' for the DuffyTools package

# one shared environment for all DuffyTools objects
DuffyToolsEnv <- new.env( hash=TRUE, parent=emptyenv())


# global constants
FASTQ_COLUMNS <- c( "READ_ID", "READ_SEQ", "SCORE")


# mapSet environment
MapSetEnv <- new.env( hash=TRUE, parent=emptyenv())
MAPSET_NAMES <- c( "speciesID", "speciesFilePrefix", "speciesText", "seqMap", "geneMap", "exonMap", "rrnaMap")


# Targets environment
TargetEnv <- new.env( hash=TRUE, parent=emptyenv())
TARGET_NAMES <- c( "TargetID", "SpeciesSet", "PrefixSet")


# Codon environment
STOP_CODON <- "*"
UNKNOWN_CODON <- "?"
CodonEnv <- new.env( hash=TRUE, parent=emptyenv())


# Alias tools
AliasEnv <- new.env( hash=TRUE, parent=emptyenv())


# LifeCycle tools
LifeCycleEnv <- new.env( hash=TRUE, parent=emptyenv())


# Ortholog tools
OrthoEnv <- new.env( hash=TRUE, parent=emptyenv())


# TimeHour tools
TimeHourEnv <- new.env( hash=TRUE, parent=emptyenv())


# Base Depth constant
EMPTY_BASE_DEPTH_TABLE <- data.frame( "START"=vector( mode="numeric", length=0), 
					"STOP"=vector( mode="numeric", length=0), 
					"DEPTH"=vector( mode="numeric", length=0))

EMPTY_BASE_DEPTH_VECTOR <- vector( mode="numeric", length=0)


# .onLoad() function is called when the package get loaded at runtime
# any needed setup goes here...
`.onLoad` <- function( libname, pkgname) {
}

`.onUnload` <- function( libpath) {
}


# .onAttach() function is called when the package gets attached,
# i.e. at the time the user first has access to the package
`.onAttach`  <- function( libname, pkgname) {

	# wake-up message
	cat( "\nPackage: \t\t", pkgname)

	# save the library and package name...
	assign( "LibraryName", value=libname, envir=DuffyToolsEnv)
	assign( "PackageName", value=pkgname, envir=DuffyToolsEnv)

	# initialize...
	DuffyTools.defaults()
}


`DuffyTools.defaults` <- function() {

	#mapset.defaults()
	target.defaults()
	codon.defaults()
}

