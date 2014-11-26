# colorTools.R  - add flexibility to R colors


`adjustColor` <- function( col, adjust=0) {
	x <- col2rgb(col); 
	dx <- if ( adjust > 0) (255 - x) else x;
	newx <- x + (dx * adjust)
	newx <- ifelse( newx > 255, 255, newx);
	newx <- ifelse( newx < 0, 0, newx);
	return( rgb( t(newx), max=255))
}


`adjustColorSet` <- function( colors, max.adjust=0.5) {

	out <- colors
	if ( any( duplicated( colors))) {
		tapply( 1:length(colors), INDEX=factor(colors), FUN=function(x) {
				if ( length(x) < 2) return()
				n <- length(x)
				delta <- max.adjust / n
				for ( i in 2:n) {
					out[x[i]] <<- adjustColor( out[x[1]], adjust=(delta*(i-1)))
				}
				return()
			})
	}
	return( out)
}


`colorBySpecies` <- function( speciesSet, palette=c('cyan', 'springgreen', 'gold', 'hotpink'),
				intergenic=NULL, intergenic.color='brown') {

	# given a vector of species IDs, return some color choices for each
	spFac <- factor( speciesSet)
	nSP <- nlevels( spFac)
	
	if ( nSP > length( palette)) palette <- rep( palette, length.out=nSP)

	colorSet <- palette[ as.numeric( spFac)]
	legendColors <- palette[ 1:nSP]
	names(legendColors) <- levels(spFac)

	if ( ! is.null( intergenic)) {
		if ( is.logical( intergenic)) intergenic <- which( intergenic)
		colorSet[ intergenic] <- intergenic.color
		legendColors <- c( legendColors, "intergenic"=intergenic.color)
	}

	out <- list( 'colors'=colorSet, 'legend.colors'=legendColors)
	return(out)
}
