# mathTools.R	-- assorted math functions


# simple moving average
# get an 'average' value for each overlapping subset of 'window' values
movingAverage <- function( x, window=9, at=names(x), by=1, FUN=mean.default, do.ends=TRUE) {

	# apply the smoothing function to a window of K adjacent points
	K <- floor((window-1)/2) * 2 + 1
	N <- length(x)

	froms <- 1:(N-K+1)
	tos <- froms + K - 1

	if ( do.ends && window > 8) {
		tosE1 <- ceiling(K/2) : (tos[1]-1)
		fromsE1 <- rep( 1, times=length(tosE1))
		fromsE2 <- ( 1:floor(K/2)) + froms[length(froms)] 
		tosE2 <- rep( N, times=length( fromsE2))
		froms <- c( fromsE1, froms, fromsE2)
		tos <- c( tosE1, tos, tosE2)
	}

	if ( by > 0) {
		use <- seq( 1, length(froms), by=round(by))
		froms <- froms[ use]
		tos <- tos[use]
	}

	values <- base::mapply( froms, tos, MoreArgs=list( "x"=x, "FUN"=FUN), FUN=function( i,j,x,FUN) {
			return( FUN( x[i:j]))
		})

	ats <- base::mapply( froms, tos, MoreArgs=list( "x"=as.numeric(at), "FUN"=mean.default), FUN=function( i,j,x,FUN) {
			return( FUN( x[i:j]))
		})

	names( values) <- ats
	return( values)
}


# takes the result of 'movingAverage' and forces the elements to have uniform distance between points
unitIntervalMA <- function( x, step=1) {

	# get the locations from the names
	atIn <- as.numeric( names(x))
	yIn <- x

	# get the integer range of X
	atRange <- range( round( atIn))
	atOut <- seq.int( atRange[1], atRange[2], by=step)

	ans <-spline( x=atIn, y=yIn, xout=atOut)

	xOut <- ans$x
	yOut <- ans$y
	# note that the spline at the very end points is a bit suspect
	yOut[1] <- yIn[1]
	names( yOut) <- xOut

	if ( xOut[ length(xOut)] < atRange[2]) {
		N <- length( yOut) + 1
		yOut[N] <- yIn[ length(yIn)]
		names(yOut)[N] <- atRange[2]
	}

	return( yOut)
}
	


# mean of a set of 'not normally distributed' values, like microarray probe intensities,
# by using a log transform, then mean, then unlog the result...
logmean <- function( x, na.rm=FALSE) {

	# zeros and negative break log
	# trap at smallest positive seen
	useX <- x
	useX[ x < 0] <- NA
	if ( sum( usable <- (! is.na(useX))) < 1) return( NA)
	if ( all( useX[ usable] == 0)) return( 0)
	isZero <- which( useX == 0)
	if ( length( isZero)) useX[ isZero] <- min( c( useX[ -isZero], 1), na.rm=T)
	logx <- log2( useX)
	m <- mean.default( logx, na.rm=na.rm)
	return( 2 ^ m)
}


logMeanPlusSD <- function( x, nSD=1, na.rm=FALSE) {

	x <- x[ x > 0]
	if ( length(x) < 1) return(0)
	logx <- log2( x)
	m <- mean.default( logx, na.rm=na.rm)
	mysd <- sd( logx, na.rm=na.rm)
	return( 2 ^ ( m + ( nSD * mysd)))
}


logMedianPlusSD <- function( x, nSD=1, na.rm=FALSE) {

	x <- x[ x > 0]
	if ( length(x) < 1) return(0)
	logx <- log2( x)
	m <- median( logx, na.rm=na.rm)
	mysd <- sd( logx, na.rm=na.rm)
	return( 2 ^ ( m + ( nSD * mysd)))
}


# square of the mean of the sqrt of X,  (i.e. generalized mean with p = 0.5)
sqrtmean <- function( x, na.rm=FALSE) {

	m <- mean.default( sqrt(x), na.rm=na.rm)
	return(m * m)
}


## Tukey biweight described here: http://www.affymetrix.com/support/technical/whitepapers/sadd_whitepaper.pdf
## "Statistical Algorithms Description Document"

# Input: a vector of values
# Output: a vector of weights based on Tukey biweight algorithm in the same order as the input

tukey.biweight = function(x, c=5, e=0.001, na.rm=FALSE){

	# c -  Tuning constant (fudge factor?)
	# e -  Small factor to prevent division by zero

	n = length( x)
	M = median( x, na.rm=na.rm)
	dM = x - M
	S = median( abs( dM), na.rm=na.rm)
	tuneConstant = rep( (c*S+e), times=n)
	u = dM / tuneConstant
	w = rep( 0, times=n)
	wTerm = (1.0 - (u * u))
	wTerm = (wTerm * wTerm)
	idx = which( abs(u)<=1);
	w[ idx] = wTerm[ idx]
	return(w);
}


tukey.mean <- function( x,  c=5,  e=0.001, na.rm=FALSE, plot=FALSE) {


	# allow a bit wider window when the number of observations is
	# too small for median to be robust...
	if (length( x) < 5) c <- c * 2

	wts <- tukey.biweight( x, c, e, na.rm=na.rm)
	m <- weighted.mean( x, wts, na.rm=na.rm)

	if ( ! plot) return( m)

	ord <- order(x)
	a <- hist( x, main="Tukey's biweight mean")
	xShow <- x[ord]
	yShow <- wts[ord]
	yMax <- max( a$counts)
	yTick <- yMax * 0.15
	yShow <- yShow * yMax
	lines( xShow, yShow, col=2, lwd=2)
	points( jitter(xShow), jitter( yShow), pch=1, cex=1)

	lines( c(m,m), c(0,yTick), col=4, lwd=3, lty=3)
	mreg <- mean(x, na.rm=T)
	lines( c(mreg,mreg), c(0,yTick), col=3, lwd=3, lty=3)
	pos <- c(2,4)
	if ( mreg < m) pos <- c(4,2)
	text( c( m, mreg), c(yTick,yTick), c( "Tukey Mean", "Regular Mean"), col=c(4,3), pos=pos)

	return( m)
}


tukey.logmean <- function( x,  c=5,  e=0.0001, na.rm=FALSE) {

	logx <- log( x, 2)
	logm <- tukey.mean( logx, c, e, na.rm=na.rm)
	m <- 2^logm
	return( m)
}


errorBar <- function( x, mode=c("se", "sd"), average.FUN=mean, plot=TRUE, at=1, whisker=0.2, error.col=1, error.lty=1) {

	mn <- average.FUN( x, na.rm=T)
	se <- s <- sd( x, na.rm=T)
	if ( match.arg( mode) == "se") {
		se <- s / sqrt( length(x))
	}

	if (plot) {
		ylo <- mn - se
		yhi <- mn + se
		lines( c( at,at), c( ylo, yhi), col=error.col, lty=error.lty)
		if ( whisker > 0) {
			lines( c( at-whisker,at+whisker), c( ylo, ylo), col=error.col, lty=error.lty)
			lines( c( at-whisker,at+whisker), c( yhi, yhi), col=error.col, lty=error.lty)
		}
	}

	return( invisible( se))
}
		

averageLineWithErrorBars <- function( x, y, average.FUN=mean, col=1, error.col=1, error.lty=1, whisker=0.2, ...) {

	xfac <- factor( x)
	ptrs <- tapply( x, xfac, FUN=NULL)
	avgY <- tapply( y, xfac, FUN=average.FUN, na.rm=T)

	lines( levels(xfac), avgY, col=col, ...)
	tapply( ptrs, xfac, function(x) {
			at <- as.numeric( levels(xfac)[ x[1]])
			myYs <- y[ ptrs == x[1]]
			errorBar( myYs, "se", average.FUN=average.FUN, at=at, whisker=whisker, 
					error.col=error.col, error.lty=error.lty)
		})

	return()
}
