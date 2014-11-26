# expressionCluster.R

# run a matrix of expression intensities  through a clustering tool

`expressionCluster` <- function( x, useLog=TRUE, FUN=diana) {


	require( cluster)

	if ( useLog) {
		x <- log( (x + 2), base=2)
	}

	clusterAns <- FUN( t(x), diss=F, metric="manhattan", stand=F, keep.diss=F, keep.data=F)

	return( clusterAns)
}

