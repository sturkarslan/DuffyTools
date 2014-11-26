# expressionMatrixToMvalule.R -- turn abundance into log 2 fold change


`expressionMatrixToMvalue` <- function( x, average.FUN=median, minIntensity=0, delta=1) {


	x <- as.matrix(x)

	if ( minIntensity < 0) minIntensity <- 0
	if (any( x <= minIntensity)) {
		cat( "Clipping low abundance values at: ", minIntensity)
		x[ x < minIntensity] <- minIntensity
	}

	rowAvgs <- apply( x, MARGIN=1, FUN=average.FUN, na.rm=T)
	isZero <- which(rowAvgs == 0)
	mv <- x
	lapply( 1:nrow(x), function(i) {
			if (i %in% isZero) {
				mv[ i, ] <<- 0
			} else {
				mv[ i, ] <<- log2( (x[i, ] + delta) / (rowAvgs[i] + delta))
			}
			return()
		})

	return( mv)
}

