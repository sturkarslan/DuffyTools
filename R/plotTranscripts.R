# plotTranscripts.R  -- assorted plots for comparison of transcripts and ratios


`plotRatios` <- function( file, geneColumn="GENE_ID", value1Column="RPKM_1_M", value2Column="RPKM_2_M", cex=1,
			units="RPKM", offset=1, keepIntergenics=FALSE, label="Plot", plotType=c("Scatter","MA"),
			marker.genes=NULL, marker.col=1, marker.cex=1, marker.labels=TRUE, marker.pch=1, minYmax=NULL, 
			sep="\t", sym.asp=TRUE, lab1="", lab2="", hideZero=FALSE) {

	tmp <- read.delim( file, as.is=T, sep=sep)
	cat( "\nRead file: ", file, "\nN_Genes: ", nrow(tmp))

	if ( !( all( c( geneColumn, value1Column, value2Column) %in% colnames(tmp)))) {
		cat( "\nMissing columns:  looked for: ", geneColumn, value1Column, value2Column, 
				"  Found: ", colnames(tmp))
		return()
	}

	# extract the parts we want
	genes <- tmp[[ geneColumn]]
	int1 <- as.numeric( tmp[[ value1Column]])
	int2 <- as.numeric( tmp[[ value2Column]])


	# allow the removal of non genes, etc.
	drops <- vector()
	if ( ! keepIntergenics) {
		drops <- grep( "(ng)", genes, fixed=TRUE)
		gmap <- getCurrentGeneMap()
		dropableGenes <- subset( gmap, SEQ_ID %in% c( "Pf3D7_PFC10_API_IRAB", "Pf3D7_M76611"))$GENE_ID
		drops2 <- which( genes %in% dropableGenes)
		drops <- base::sort( base::union( drops, drops2))
	}
	if ( length( drops) > 0) {
		genes <- genes[ -drops]
		int1 <- int1[ -drops]
		int2 <- int2[ -drops]
		cat( "\nAfter dropping non-genes: ", length(genes))
	}

	plotType <- base::match.arg( plotType)
	if ( plotType == "MA") {
		ans <- makeMAplot( genes, int1, int2, offset=offset, units=units, label=label, cex=cex,
			marker.genes=marker.genes, marker.col=marker.col, marker.cex=marker.cex,
			marker.labels=marker.labels, marker.pch=marker.pch, 
			minYmax=minYmax, hideZero=hideZero)
	} else {
		ans <- makeScatterplot( genes, int1, int2, offset=offset, units1=units, label=label, cex=cex,
			marker.genes=marker.genes, marker.col=marker.col, marker.cex=marker.cex,
			marker.labels=marker.labels, marker.pch=marker.pch, 
			minYmax=minYmax, sym.asp=sym.asp, lab1=lab1, lab2=lab2, hideZero=hideZero)
	}

	return( invisible( ans))
}


`plotTranscripts` <- function( file1, file2, geneColumn="GENE_ID", value1Column="RPKM_M", value2Column=value1Column, cex=1,
			units1="RPKM", units2=units1, offset=1, keepIntergenics=FALSE, label="Plot", plotType=c("Scatter","MA"),
			marker.genes=NULL, marker.col=1, marker.cex=1, marker.labels=TRUE, marker.pch=1, minYmax=NULL, 
			sep="\t", sym.asp=TRUE, lab1="", lab2="", hideZero=FALSE) {

	tmp <- read.delim( file1, as.is=T, sep=sep)
	cat( "\nRead file: ", file1, "\nN_Genes: ", nrow(tmp), "\n")
	if ( !( all( c( geneColumn, value1Column) %in% colnames(tmp)))) {
		cat( "\nMissing columns:  looked for: ", geneColumn, value1Column, "  Found: ", colnames(tmp))
		return()
	}
	# extract the parts we want
	genes1 <- tmp[[ geneColumn]]
	int1 <- as.numeric( tmp[[ value1Column]])

	tmp <- read.delim( file2, as.is=T, sep=sep)
	cat( "\nRead file: ", file2, "\nN_Genes: ", nrow(tmp), "\n")
	if ( !( all( c( geneColumn, value2Column) %in% colnames(tmp)))) {
		cat( "\nMissing columns:  looked for: ", geneColumn, value2Column, "  Found: ", colnames(tmp))
		return()
	}
	# extract the parts we want
	genes2 <- tmp[[ geneColumn]]
	int2 <- as.numeric( tmp[[ value2Column]])

	# combine and resolve
	both <- intersect( genes1, genes2)
	wh1 <- base::match( both, genes1)
	genes1 <- genes1[ wh1]
	int1 <- int1[ wh1]
	wh2 <- base::match( both, genes2)
	genes2 <- genes2[ wh2]
	int2 <- int2[ wh2]
	genes <- both

	# put into MA order
	v <- log2( (int1+1)/(int2+1))
	ord <- base::order( v, decreasing=T)
	genes <- genes[ord]
	int1 <- int1[ord]
	int2 <- int2[ord]

	# allow the removal of non genes, etc.
	drops <- vector()
	if ( ! keepIntergenics) {
		drops <- grep( "(ng)", genes, fixed=TRUE)
		gmap <- getCurrentGeneMap()
		dropableGenes <- subset( gmap, SEQ_ID %in% c( "Pf3D7_PFC10_API_IRAB", "Pf3D7_M76611"))$GENE_ID
		drops2 <- which( genes %in% dropableGenes)
		drops <- base::sort( base::union( drops, drops2))
	}
	if ( length( drops) > 0) {
		genes <- genes[ -drops]
		int1 <- int1[ -drops]
		int2 <- int2[ -drops]
		cat( "\nAfter dropping non-genes: ", length(genes))
	}
	cat( "\nN_Genes in common: ",  length(genes), "\n")

	plotType <- match.arg( plotType)
	if ( plotType == "MA") {
		ans <- makeMAplot( genes, int1, int2, offset=offset, units=units1, label=label, cex=cex,
			marker.genes=marker.genes, marker.col=marker.col, marker.cex=marker.cex,
			marker.labels=marker.labels, marker.pch=marker.pch, 
			minYmax=minYmax, hideZero=hideZero)
	} else {
		ans <- makeScatterplot( genes, int1, int2, offset=offset, units1=units1, units2=units2, 
			label=label, cex=cex,
			marker.genes=marker.genes, marker.col=marker.col, marker.cex=marker.cex,
			marker.labels=marker.labels, marker.pch=marker.pch, 
			minYmax=minYmax, sym.asp=sym.asp, lab1=lab1, lab2=lab2, hideZero=hideZero)
	}

	return( invisible( ans))
}


`plotRatiosTwoFiles` <- function( file1, file2, geneColumn="GENE_ID", valueColumn="LOG2FOLD_M", cex=1,
			units="Log2 Fold", offset=0, keepIntergenics=FALSE, label="Plot", 
			marker.genes=NULL, marker.col=1, marker.cex=1, marker.labels=TRUE, marker.pch=1, minYmax=NULL, 
			sep="\t", sym.asp=TRUE, lab1="", lab2="", hideZero=FALSE) {

	tmp <- read.delim( file1, as.is=T, sep=sep)
	cat( "\nRead file: ", file1, "\nN_Genes: ", nrow(tmp), "\n")
	if ( !( all( c( geneColumn, valueColumn) %in% colnames(tmp)))) {
		cat( "\nMissing columns:  looked for: ", geneColumn, valueColumn, "  Found: ", colnames(tmp))
		return()
	}
	# extract the parts we want
	genes1 <- tmp[[ geneColumn]]
	int1 <- as.numeric( tmp[[ valueColumn]])

	tmp <- read.delim( file2, as.is=T, sep=sep)
	cat( "\nRead file: ", file2, "\nN_Genes: ", nrow(tmp), "\n")
	if ( !( all( c( geneColumn, valueColumn) %in% colnames(tmp)))) {
		cat( "\nMissing columns:  looked for: ", geneColumn, valueColumn, "  Found: ", colnames(tmp))
		return()
	}
	# extract the parts we want
	genes2 <- tmp[[ geneColumn]]
	int2 <- as.numeric( tmp[[ valueColumn]])

	# combine and resolve
	both <- intersect( genes1, genes2)
	wh1 <- base::match( both, genes1)
	genes1 <- genes1[ wh1]
	int1 <- int1[ wh1]
	wh2 <- base::match( both, genes2)
	genes2 <- genes2[ wh2]
	int2 <- int2[ wh2]
	genes <- both

	ord <- base::order( int1, decreasing=T)
	genes <- genes[ord]
	int1 <- int1[ord]
	int2 <- int2[ord]

	# allow the removal of non genes, etc.
	drops <- vector()
	if ( ! keepIntergenics) {
		drops <- grep( "(ng)", genes, fixed=TRUE)
		gmap <- getCurrentGeneMap()
		dropableGenes <- subset( gmap, SEQ_ID %in% c( "Pf3D7_PFC10_API_IRAB", "Pf3D7_M76611"))$GENE_ID
		drops2 <- which( genes %in% dropableGenes)
		drops <- base::sort( base::union( drops, drops2))
	}
	if ( length( drops) > 0) {
		genes <- genes[ -drops]
		int1 <- int1[ -drops]
		int2 <- int2[ -drops]
		cat( "\nAfter dropping non-genes: ", length(genes))
	}
	cat( "\nN_Genes in common: ",  length(genes), "\n")

	ans <- makeScatterplot( genes, int1, int2, offset=offset, units1=units, label=label, cex=cex,
		marker.genes=marker.genes, marker.col=marker.col, marker.cex=marker.cex,
		marker.labels=marker.labels, marker.pch=marker.pch, 
		minYmax=minYmax, sym.asp=sym.asp, lab1=lab1, lab2=lab2, useLog="", hideZero=hideZero)
	
	return( invisible( ans))
}


`makeMAplot` <- function( genes, int1, int2, offset=1, units="Intensity", label="MA plot", 
			marker.genes=NULL, marker.col=1, marker.cex=1, marker.labels=TRUE, 
			marker.pch=1, minYmax=NULL, cex=1, hideZero=FALSE) {

	if (hideZero) {
		zeros <- which( int1 == 0 | int2 == 0)
		if( length(zeros) > 0) {
			int1 <- int1[ -zeros]
			int2 <- int2[ -zeros]
			genes <- genes[ -zeros]
		}
	}

	# allow a linear shift of background intensity
	int1 <- int1 + offset
	int2 <- int2 + offset

	# A is average, M is fold change
	a <- log( sqrt(int1 * int2), 2)
	m <- log( (int1 / int2), 2)
	myYrange <- range(m) * 1.2
	if ( !is.null(minYmax)) {
		if (myYrange[2] < minYmax) myYrange[2] <- minYmax
		if (myYrange[1] > (-minYmax)) myYrange[1] <- (-minYmax)
	}

	plot ( a, m, main=label, xlab=paste( "A  log_2( Average",units,")"), 
		ylab=paste("M  log_2( Ratio",units,")"), ylim=myYrange, cex=cex, font.axis=2, font.lab=2)

#	addConfidenceLines( a, m, conf=0.95)

	# label selected genes, assumes the genes 'up in set1' are in the first half...
	if ( length( marker.genes) > 0) {
		if ( marker.genes[1] == "identify") {
			identify( a, m, shortGeneName(genes, keep=1), col=marker.col, cex=marker.cex)
		} else {
		who <- base::match( marker.genes, genes, nomatch=0)
		marker.genes <- marker.genes[ who > 0]
		if ( length( marker.col) > 1) marker.col <- marker.col[ who > 0]
		who <- who[ who > 0]
		if ( length(who) > 0) {
			# put name above for first half, below for second half...
			pos <- ifelse( m[who] > 0, 3, 1)
			points( a[who], m[who], col=marker.col, bg=marker.col, pch=marker.pch, cex=cex)
			if (marker.labels) text( a[who], m[who], genes[who], pos=pos, col=marker.col, cex=marker.cex)
		}}
	}

	Rpearson <- cor( int1, int2, use="complete", method="pearson")
	Rspearman <- cor( int1, int2, use="complete", method="spearman")

	return( list( "x"=a, "y"=m, "id"=genes, "Pearson_R"=Rpearson, "Spearman_Rho"=Rspearman))
}


`makeScatterplot` <- function( genes, int1, int2, offset=1, units1="Intensity", 
			units2=units1, label="MA plot", 
			marker.genes=NULL, marker.col=1, marker.cex=1, marker.labels=TRUE, 
			marker.pch=1, minYmax=NULL, lab1="", lab2="", cex=1, useLog="xy", 
			hideZero=FALSE, sym.asp=TRUE, show.cor=TRUE) {

	# we are plotting file 1 on Y axis, file 2 on X axis

	if (hideZero) {
		zeros <- which( int1 == 0 | int2 == 0)
		if( length(zeros) > 0) {
			int1 <- int1[ -zeros]
			int2 <- int2[ -zeros]
			genes <- genes[ -zeros]
		}
	}

	# allow a linear shift of background intensity
	int1 <- int1 + offset
	int2 <- int2 + offset

	myRange1 <- range( int1)
	myRange2 <- range( int2)
	if ( !is.null(minYmax)) {
		if (myRange1[2] < minYmax) myRange1[2] <- minYmax
		if (myRange2[2] < minYmax) myRange2[2] <- minYmax
	}
	if (sym.asp) myRange1 <- myRange2 <- range( c( myRange1, myRange2))

	plot ( int2, int1, main=label, log=useLog, xlab=paste( lab2, "    (",units2, ")", sep=""), 
		ylab=paste( lab1, "    (", units1, ")", sep=""), xlim=myRange2, ylim=myRange1, 
		cex=cex, font.axis=2, font.lab=2, xaxt="n", yaxt="n")

	xAts <- c( 0, axTicks(1)+offset)
	xAtLbls <- c( 0, axTicks(1))
	axis( 1, at=xAts, label=xAtLbls, font=2)
	yAts <- c( 0, axTicks(2)+offset)
	yAtLbls <- c( 0, axTicks(2))
	axis( 2, at=yAts, label=yAtLbls, font=2)

	if (show.cor) {
		cat( "\nPearson's R =", cor( int1, int2), "\nN_Zeros: ", sum(int1 <= offset), sum(int2 <= offset))

		if (FALSE) {
		toUse <- which(  int1 > offset & int2 > offset)
		Xsave <- log2(int2[toUse])
		Ysave <- log2(int1[toUse])
		lsAns <- lsfit( x=Xsave, y=Ysave)
		orig <- coef(lsAns)[1]
		slope <- coef(lsAns)[2]
		mySteps <- c( -1000, -100, -10, -1, 0, 1, 10, 100, 1000, 10000, 100000)
		myX <- mySteps
		myY <- orig + (mySteps * slope)
		lines( (2^myX), (2^myY), col=2, lwd=2)
		#addConfidenceLines( int2, int1, conf=0.95)
		}
	}

	# label selected genes, assumes the genes 'up in set1' are in the first half...
	if ( length( marker.genes) > 0) {
		if ( marker.genes[1] == "identify") {
			identify( int2, int1, shortGeneName(genes, keep=1), col=marker.col, cex=marker.cex)
		} else {
		who <- base::match( marker.genes, genes, nomatch=0)
		marker.genes <- marker.genes[ who > 0]
		if ( length( marker.col) > 1) marker.col <- marker.col[ who > 0]
		who <- who[ who > 0]
		if ( length(who) > 0) {
			# put name left for first half, right for second half...
			pos <- ifelse( int1[who] > int2[who], 2, 4)
			points( int2[who], int1[who], col=marker.col, bg=marker.col, pch=marker.pch, cex=cex)
			if (marker.labels) text( int2[who], int1[who], genes[who], pos=pos, col=marker.col, cex=marker.cex)
		}}
	}

	Rpearson <- cor( int1, int2, use="complete", method="pearson")
	Rspearman <- cor( int1, int2, use="complete", method="spearman")
	return( list( "x"=int2, "y"=int1, "id"=genes, "Pearson_R"=Rpearson, "Spearman_Rho"=Rspearman))
}


addConfidenceLines <- function( int1, int2, conf=0.95) {

	# get the regression line model, and its fitted values with std errs.
	model <- lm( int2 ~ int1)
	preds <- predict( model, se.fit=TRUE)

	# 'Working-Hotelling multiplier'
	W <- sqrt( 2 * qf( conf, 2, (length(int1) - 2)))

	# upper and lower conf bands
	dw <- W * preds$se.fit
	bands <- data.frame( "lowBand"=preds$fit - dw, "highBand"=preds$fit + dw)

	# the regression line
	abline( model$coeficents[1], model$coeficents[2], untf=TRUE, col=2, lty=1)

	ordX <- base::order(int1)
	lines( int1[ordX], bands$lowBand[ordX], lty=2, col=2)
	lines( int1[ordX], bands$highBand[ordX], lty=2, col=2)
	return()
}
