# pipe.RadarPlots.R -- visualize gene DE as a family of radar (spider) plots

`pipe.RadarPlots` <- function( sampleIDset, speciesID="Pf3D7", annotationFile="Annotation.txt",
				optionsFile="Options.txt", results.path=NULL, folderName=NULL,
				tool=c( "MetaResults", "DESeq", "EdgeR", "RankProduct", "Roundrobin", "SAM"),
				groupColumn="Group", colorColumn="Color", 
				geneSets=c("GO.BiologicalProcess", "GO.MolecularFunction",
				"GO.CellularComponent", "KEGG.Pathways", "MetabolicPathways", 
				"GeneProduct", "PBMC.GeneModules", "Blood.GeneModules"), 
				legend.prefix=NULL, legend.order=NULL, Nshow=36, 
				start=pi/4, radial.labels=FALSE, radial.margin=c( 2,2,6,2),
				radial.lim=NULL, boxed.radial=F, label.prop=1, lwd=5, main=NULL, ...)
{

	require( plotrix)
	if ( names(dev.cur()) != "X11") x11( bg='white', width=12, height=8)

	setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	# get the SampleIDs that we know are valid
	annT <- readAnnotationTable( annotationFile)
	allSamples <- sampleIDset
	myAnnT <- subset( annT, SampleID %in% allSamples)
	allSamples <- myAnnT$SampleID
	if ( ! (groupColumn %in% colnames(myAnnT))) stop( paste( "Given grouping column not in annotation table: ", groupColumn))
	if ( ! (colorColumn %in% colnames(myAnnT))) stop( paste( "Given coloring column not in annotation table: ", colorColumn))

	# build the path to the DE results we will use, and put our Radar results
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", verbose=F)
	}
	if ( is.null( folderName)) stop( "Explicit Folder Name of DE results is required")
	folderName <- paste( prefix, folderName, sep=".")
	tool <- match.arg( tool)
	dePath <- file.path( results.path, tool, folderName)
	if ( ! file.exists( dePath)) stop( paste( "DE results folder not found: ", dePath))
	radarPath <- file.path( dePath, "RadarPlots")
	if ( ! file.exists( radarPath)) dir.create( radarPath, recursive=T)

	# we have to read the transcriptomes in... to pre-load the data
	needLoad <- TRUE
	if (needLoad) {
		cat( "\nLoading ", length( allSamples), "transcriptomes..")
		files <- paste( allSamples, prefix, "Transcript.txt", sep=".")
		files <- file.path( results.path, "transcript", files)
		radarM <<- expressionFileSetToMatrix( files, allSamples, verbose=T)
		cat( "\nConverting Expression Abundance to M-values..")
		radarMA <<- expressionMatrixToMvalue( radarM)
		cat( "  Done.\n")
	}


	# now we are read to make those radar plots
	if ( is.list( geneSets)) {
		ans <- pipe.OneRadarPlot( sampleIDset=allSamples, speciesID=speciesID, annotationFile=annotationFile,
				optionsFile=optionsFile, results.path=results.path, groupColumn=groupColumn,
				colorColumn=colorColumn, geneSet=geneSets, reload=FALSE, Nshow=Nshow, start=start, 
				radial.labels=radial.labels, radial.margin=radial.margin, radial.lim=radial.lim,
				boxed.radial=boxed.radial, label.prop=label.prop, lwd=lwd, main=main, ...)
		geneSetName <- "Radar"
		plotFile <- file.path( radarPath, paste( geneSetName, "png", sep="."))
		dev.print( png, plotFile, width=900, height=640)
		csvFile <- file.path( radarPath, paste( geneSetName, "csv", sep="."))
		write.table( ans, csvFile, sep=",", quote=T, row.names=F)

	} else {
		for (gs in geneSets) {
			ans <- pipe.OneRadarPlot( sampleIDset=allSamples, speciesID=speciesID, annotationFile=annotationFile,
					optionsFile=optionsFile, results.path=results.path, groupColumn=groupColumn,
					colorColumn=colorColumn, geneSet=gs, reload=FALSE, Nshow=Nshow, start=start, 
					radial.labels=radial.labels, radial.margin=radial.margin, radial.lim=radial.lim,
					boxed.radial=boxed.radial, label.prop=label.prop, lwd=lwd, main=main, ...)

			plotFile <- file.path( radarPath, paste( "Radar", gs, "png", sep="."))
			dev.print( png, plotFile, width=900, height=640)
			csvFile <- file.path( radarPath, paste( "Radar", gs, "csv", sep="."))
			write.table( ans, csvFile, sep=",", quote=T, row.names=F)
		}
	}

	return()
}


`pipe.OneRadarPlot` <- function( sampleIDset, speciesID="Pf3D7", annotationFile="Annotation.txt",
				optionsFile="Options.txt", results.path=NULL,  
				groupColumn="Group", colorColumn="Color", 
				legend.prefix=NULL, legend.order=NULL,
				geneSet=c("GO.BiologicalProcess", "GO.MolecularFunction",
				"GO.CellularComponent", "KEGG.Pathways", "MetabolicPathways", 
				"GeneProduct", "PBMC.GeneModules", "Blood.GeneModules", "GSEA"), 
				reload=FALSE, Nshow=36, 
				start=pi/4, radial.labels=FALSE, radial.margin=c( 2,2,6,2),
				radial.lim=NULL, boxed.radial=F, label.prop=1, lwd=5, main=NULL, ...)
{

	setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	mainText <- "Radar Plot:  "
	if ( is.character( geneSet)) {
		if ( missing(geneSet)) geneSet <- match.arg( geneSet)
		geneSetFile <- paste( prefix, geneSet, sep=".")
		data( list=geneSetFile)
		if ( ! exists( "allGeneSets")) stop( paste( "Failed to load dataset: ", geneSetFile))
		mainText <- paste( mainText, "    ", geneSet)
		if ( ! is.null( main)) mainText <- paste( mainText, main, sep="\n")
	} else {
		if ( typeof( geneSet) != "list") stop( "'geneSet' must be a character string or a list")
		allGeneSets <- geneSet
		if ( ! is.null( main)) mainText <- paste( mainText, main, sep="    ")
	}

	# use the samples and grouping column to know the colors and files to load
	annT <- readAnnotationTable( annotationFile)
	allSamples <- sampleIDset
	myAnnT <- subset( annT, SampleID %in% allSamples)
	allSamples <- myAnnT$SampleID

	if ( ! (groupColumn %in% colnames(myAnnT))) stop( paste( "Given grouping column not in annotation table: ", groupColumn))
	grpFac <- factor( myAnnT[[ groupColumn]])
	allGroups <- tapply( myAnnT$SampleID, grpFac, FUN=c)
	NperGroup <- sapply( allGroups, length)

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", verbose=F)
	}

	# map from the group names to find the colors to use...
	if ( ! (colorColumn %in% colnames(myAnnT))) stop( paste( "Given coloring column not in annotation table: ", colorColumn))
	where <- base::match( names(allGroups), myAnnT[[ groupColumn]])
	mycolors <- myAnnT[[ colorColumn]][ where]

	# see if we have to read the transcriptomes in...
	needLoad <- TRUE
	if ( exists( "radarM") && all( colnames(radarM) == allSamples) && nrow(radarM) > 100 && !reload) needLoad <- FALSE
	if (needLoad) {
		files <- paste( allSamples, prefix, "Transcript.txt", sep=".")
		files <- file.path( results.path, "transcript", files)
		radarM <<- expressionFileSetToMatrix( files, allSamples, verbose=T)
		radarMA <<- expressionMatrixToMvalue( radarM)
	}

	# try to standarize the names of the gene sets to keep them easy to view
	fullPathNames <- names(allGeneSets)
	names(allGeneSets) <- cleanModuleNames( names(allGeneSets))

	# do the reduction & grouping
	radarAns <- reduceMatrixToModules( radarMA, geneModules=allGeneSets, sampleTraits=allGroups,
				gene.names=shortGeneName( rownames( radarMA), keep=1), sample.names=colnames(radarMA))
	mShow <- radarMOD <- radarAns$matrix
	pShow <- radarPvalue <- radarAns$p.value

	# trim down to less groups if we need to
	Nshow <- min( Nshow, nrow( radarMOD))

	# use both Pvalue and magnitudes
	magnitudes <- diff( apply( radarMOD, 1, range))
	bestPs <- apply( radarPvalue, 1, min)
	ord <- diffExpressRankOrder( magnitudes, bestPs)
	mShow <- radarMOD[ ord[1:Nshow], ]
	pShow <- radarPvalue[ ord[1:Nshow], ]
	mOut <- radarMOD[ ord, ]
	pOut <- radarPvalue[ ord, ]
	fullPathNames <- fullPathNames[ ord]

	# now alphabetize on any part after the Module names
	ord <- order( sub( "M.+: +", "", rownames(mShow)))
	mShow <- mShow[ ord, ]
	pShow <- pShow[ ord, ]

	# plot it now
	require( plotrix)
	if ( is.null( radial.lim)) radial.lim <- range( as.vector( mShow)) * 1.5
	radial.plot( t(mShow), labels=rownames(mShow), radlab=radial.labels, rp.type="p", line.col=mycolors,
			start=start, clockwise=T, mar=radial.margin, radial.lim=radial.lim, label.prop=label.prop,
			show.grid.labels=3, lwd=lwd, main=mainText, ...)

	# take more control of the legend location
	usr <- par( "usr")

	legendText <- names(allGroups)
	if ( any( NperGroup > 1)) legendText <- paste( names(allGroups), "  (N=",NperGroup, ")", sep="")
	if ( ! is.null(legend.prefix)) legendText <- paste( legend.prefix, legendText, sep=":  ")
	if ( ! is.null(legend.order)) {
		legendText <- legendText[ legend.order]
		mycolors <- mycolors[ legend.order]
	}
	legend( x=usr[1]*1.5, y=usr[4], legendText, lwd=lwd, col=mycolors, bg="white", cex=1.1)

	# return the full table, what we drew is at the very top
	out <- data.frame( "PathName"=fullPathNames, round(mOut,digits=4), formatC( pOut,format="e",digits=2),
				stringsAsFactors=F)
	colnames(out) <- c( "PathName", paste( "Fold", colnames(mOut), sep="_"), paste( "Pvalue", colnames(mOut), sep="_"))
	rownames(out) <- 1:nrow(out)

	return(out)
}


`cleanModuleNames` <- function( nams) {

	out <- nams

	# 1. strip out any hyperlink anchors
	out <- sub( " ?<a.+/a>", "", out)

	# 2.  any () get a second line and may get clipped
	hasParen <- grep( "(", out, fixed=T)
	fronts <- sub( "\\(.+", "", out[hasParen])
	backs <- sub( "(.+)(\\(.+)", "\\2", out[hasParen])
	doClip <- which( nchar( backs) > 40)
	backs[ doClip] <- paste( substr( backs[doClip], 1, 37), "...)", sep="")
	out[ hasParen] <- paste( fronts, backs, sep="\n")

	out
}



`radarPlot` <- function( m, row.names=rownames(m), colGroups=colnames(m), colors=1:ncol(m),
				legend.prefix=NULL, legend.order=NULL, foldChangeTransform=TRUE,
				geneSet=c("GO.BiologicalProcess", "GO.MolecularFunction",
				"GO.CellularComponent", "KEGG.Pathways", "MetabolicPathways", 
				"GeneProduct", "PBMC.GeneModules", "Blood.GeneModules", "GSEA"), 
				reload=FALSE, Nshow=36, 
				start=pi/4, radial.labels=FALSE, radial.margin=c( 2,2,6,2),
				radial.lim=NULL, boxed.radial=F, label.prop=1, lwd=5, main=NULL, ...)
{

	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.character( geneSet)) {
		geneSet <- match.arg( geneSet)
		geneSetFile <- paste( prefix, geneSet, sep=".")
		data( list=geneSetFile)
		if ( ! exists( "allGeneSets")) stop( paste( "Failed to load dataset: ", geneSetFile))
		if (is.null( main)) main <- paste( "Radar Plot:  ", geneSet)
	} else {
		if ( typeof( geneSet) != "list") stop( "'geneSet' must be a character string or a list")
		allGeneSets <- geneSet
		if (is.null( main)) main <- "Radar Plot:  "
	}

	allGroups <- colGroups
	NperGroup <- sapply( allGroups, length)
	mycolors <- colors[1:length(allGroups)]

	if ( foldChangeTransform) {
		radarM <<- m
		radarMA <<- expressionMatrixToMvalue( radarM)
	} else {
		radarMA <<- m
	}

	# try to standarize the names of the gene sets to keep them easy to view
	fullPathNames <- names(allGeneSets)
	names(allGeneSets) <- cleanModuleNames( names(allGeneSets))

	# do the reduction & grouping
	radarAns <- reduceMatrixToModules( radarMA, geneGroups=allGeneSets, sampleGroups=allGroups,
						row.names=row.names)
	mShow <- radarMOD <- radarAns$matrix
	pShow <- radarPvalue <- radarAns$p.value

	# trim down to less groups if we need to
	Nshow <- min( Nshow, nrow( radarMOD))

	# use both Pvalue and magnitudes
	magnitudes <- diff( apply( radarMOD, 1, range))
	bestPs <- apply( radarPvalue, 1, min)
	ord <- diffExpressRankOrder( magnitudes, bestPs)
	mShow <- radarMOD[ ord[1:Nshow], ]
	pShow <- radarPvalue[ ord[1:Nshow], ]
	mOut <- radarMOD[ ord, ]
	pOut <- radarPvalue[ ord, ]
	fullPathNames <- fullPathNames[ ord]

	# alphabetize on any part after the Module names
	ord <- order( sub( "M.+: +", "", rownames(mShow)))
	mShow <- mShow[ ord, ]
	pShow <- pShow[ ord, ]

	# plot it now
	require( plotrix)
	if ( is.null( radial.lim)) radial.lim <- range( as.vector( mShow)) * 1.5
	radial.plot( t(mShow), labels=rownames(mShow), radlab=radial.labels, rp.type="p", line.col=mycolors,
			start=start, clockwise=T, mar=radial.margin, radial.lim=radial.lim, label.prop=label.prop,
			show.grid.labels=3, lwd=lwd, main=main, ...)

	# take more control of the legend location
	usr <- par( "usr")

	legendText <- paste( names(allGroups), "  (N=",NperGroup, ")", sep="")
	if ( ! is.null(legend.prefix)) legendText <- paste( legend.prefix, legendText, sep=":  ")
	if ( ! is.null(legend.order)) {
		legendText <- legendText[ legend.order]
		mycolors <- mycolors[ legend.order]
	}
	legend( x=usr[1]*1.35, y=usr[4], legendText, lwd=lwd, col=mycolors, bg="white", cex=1.1)

	# return what we drew
	out <- data.frame( "PathName"=fullPathNames, mOut, pOut, stringsAsFactors=F)
	colnames(out) <- c( "PathName", paste( "Fold", colnames(mOut), sep="_"), paste( "Pvalue", colnames(mOut), sep="_"))
	rownames(out) <- 1:nrow(out)

	radarAns <<- out

	return(out)
}
