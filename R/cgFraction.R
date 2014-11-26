`cgFraction` <- function(dna) {

	# get the first element of the list returned from "strsplit"
	letters <- strsplit( toupper( dna[1]), split="", fixed=TRUE)[[1]]
	#n <- length(letters)
	#return( as.double(c) / as.double(n) )
	at <- sum( letters %in% c( "A","T"))
	cg <- sum( letters %in% c( "C","G"))
	return( as.double(cg) / as.double(at+cg))
}

