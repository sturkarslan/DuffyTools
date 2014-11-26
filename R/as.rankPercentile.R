# as.rankPercentile.R

`as.rankPercentile` <- function( x) {

	N <- length(x)
	ord <- order( x, na.last=FALSE)
	ranks <- rank( x, na.last=T, ties.method="average")
	rankPct <- vector( mode="numeric", length=N)
	#rankPct[ord] <- (1:N * 100 / N)
	rankPct <- ( ranks * 100 / N)
	return( rankPct)
}
