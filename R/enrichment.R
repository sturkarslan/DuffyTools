# enrichment.R

# calcualate the probability of seeing a certain number of genes in common
`enrichment` <- function( nMatch, nYourSet, nTotal, nTargetSubset) {

	x <- nMatch
	m <- nTargetSubset
	n <- nTotal - m
	k <- nYourSet

	expected <- (nTargetSubset / nTotal) * nYourSet

	# get the entire prob. dist. and sum up both halves
	allPs <- dhyper( 0:k, m, n, k)
	# the probability of 0 to X genes
	lowerTailPvalue <- sum( allPs[1:(x+1)])
	# the probability of X up to K genes
	upperTailPvalue <- sum( allPs[(x+1):length(allPs)])

	out <- list( "nWhiteBalls"=m, "nBlackBalls"=n, "nDrawn"=k, "nWhiteDrawn"=x, 
			"nExpected"=expected, "P_atLeast_N"=upperTailPvalue, "P_atMost_N"=lowerTailPvalue)

	return( out)
}


enrichment.Nway <- function( nSets=2, nMatch, nDrawn, nTotal, nSimulations=1000000) {

	nTotal <- as.integer( nTotal)
	nDrawn <- as.integer( nDrawn)

	tallyEnrichment <- function( nOverlap, nGood, nSimulation) {

		# trim storage and tabulate
		if ( length( nOverlap) > nGood) length( nOverlap) <- nGood
		dist <- table( nOverlap)
		expected <- sum( nOverlap) / nSimulation

		yes <-sum( nOverlap >= nMatch)
		no <- nSimulations - yes

		return( list( "nExpected"=expected, 
			"P_atLeast_N"=( yes / nSimulations), 
			"P_atMost_N"=( no / nSimulations), 
			"distribution"=dist))
	}

	# storage to hold how many in common from all successful trials
	nOverlap <- vector( mode="numeric", length=nSimulations/10)
	nGood <- 0

	# do those trials
	base::lapply( 1:nSimulations, function(x) {
		picks <- base::unlist( base::lapply( 1:nSets, FUN=function(x) sample.int( n=nTotal, size=nDrawn)))
		hits <- sum( tabulate(picks) == nSets)
		if ( hits > 0) {
			nGood <<- nGood + 1
			nOverlap[ nGood] <<- hits
		}

		if ( x %% 10000 == 0) {
			ans <- tallyEnrichment( nOverlap, nGood, x)
			cat( "\rIter:", x, "  Expect:", ans$nExpected, "  P.atleastN:", ans$P_atLeast_N,
					"  P.atmostN:", ans$P_atMost_N,
					"  Dist: ", paste( names( ans$distribution), ans$distribution, sep="-",
					collapse=", "))
		}
		return()
	})

	# trim storage and tabulate
	return( tallyEnrichment( nOverlap, nGood, nSimulations))
}


simulate.enrichment.Nway <- function( nSets=3, nMatch=40, nDrawn=100, nTotal=5475, nSimulations=500000) {

	cat( "\nChoose", nDrawn, "from a pool of", nTotal, "into", nSets, "sets.")
	cat( "\nLikelihood of finding ", nMatch, "in common:")
	cat( "\nSimulating ", nSimulations, " random trials...")
	ans <- enrichment.Nway( nSets=nSets, nMatch=nMatch, nDrawn=nDrawn, nTotal=nTotal, nSimulations=nSimulations)
	cat( "\n")
	print( ans)

	hits <<- ans$distribution
	N <<- length( hits)
	pcts <- vector()
	for  ( i in 1:N)  pcts[i] <- sum( hits[i:N]) / nSimulations
	log10pcts <- log10( pcts)
	
	fitans <<- lsfit( 0:N, c( 0, log10pcts))

	plot( 0:N, c(0,log10pcts), xlim=c(0,max(20, N)), ylim=c(-20,0), 
		main=paste( "Probability Curve of 'Choose ", nDrawn, " from ",nTotal,"' overlap in ",nSets," datasets"), 
		xlab=paste("Number of genes found in 'Top ",nDrawn,"' of all ",nSets," datasets"), ylab= "log_10( P )")

	abline( fitans, col=2)
	dev.print( png, "./test.png", width=700, height=500)
	return( ans)
}


estimate.enrichment.Nway.Pvalue <- function( nMatch=5) {

	yEstimate <- coef(fitans)[1] + (coef(fitans)[2] * nMatch)
	return( as.vector( 10 ^ yEstimate))
}
