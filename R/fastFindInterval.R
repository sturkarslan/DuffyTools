# fastFindInterval.R

#  bare bones version of 'findInterval()' for speed
`fastFindInterval` <- function( x, vec) {

	return( findInterval( x, vec))

	nx <- length(x)
	index <- integer(nx)
	.C( "find_interv_vec", nt=as.double(vec), n=length(vec), x=as.double(x), nx=nx, FALSE, FALSE,
		index, DUP=FALSE, NAOK=TRUE, PACKAGE="base")
	return( index)
}
