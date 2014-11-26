# tmTools.R

`calcTM` <- function( probes) {

	# put a Perl script in the current folder...
	perlFile <- "./MeltingTemp.pl"
	if ( ! file.exists( perlFile)) {
		data( MeltingTempScript, envir=environment())
		writeLines( perlScript, perlFile)
		system( paste( "chmod +x ", perlFile))
	}

	probeSet <- toupper( probes)

	tmSet <- as.numeric( system( perlFile, intern=TRUE, input=probeSet))

	return( tmSet)
}

