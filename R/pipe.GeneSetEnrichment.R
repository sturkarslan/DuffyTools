# pipe.GeneSetEnrichment.R -- investigate sets of DE genes like pathway and GO groups by enrichment

`pipe.GeneSetEnrichment` <- function( sampleIDset, speciesID="Pf3D7", annotationFile="Annotation.txt",
				optionsFile="Options.txt", results.path=NULL,  folderName="", 
				toolName=c( "MetaResults", "RoundRobin", "RankProduct", "SAM", "EdgeR"), 
				geneColumn=if(speciesID %in% c("Hs_grc","MacMu","Mmu_grc")) "GENE_NAME" else "GENE_ID", 
				groupColumn="Group", 
				geneSets=c("GO.BiologicalProcess", "GO.MolecularFunction",
				"GO.CellularComponent", "KEGG.Pathways", "MetabolicPathways", 
				"GeneProduct", "PBMC.GeneModules", "Blood.GeneModules"), 
				descriptor="Enrichment", maxPvalue=0.05, wt.enrich=1, wt.pvalue=2,
				verbose=TRUE) {

	toolName <- match.arg( toolName)

	if (verbose) {
		cat( "\n\n=============================")
		cat( "\nStarting pipe 'GeneSetEnrichment' on Sample Set: \n")
		print(sampleIDset)
		cat(  "\nUsing results from Species:  ", speciesID)
		cat(  "\nUsing results from DE tool:  ", toolName)
		cat(  "\nGene Sets to analyze:        ", geneSets,"\n")
	}

	annT <- readAnnotationTable( annotationFile)
	allSamples <- unique( unlist( sampleIDset))
	myAnnT <- subset( annT, SampleID %in% allSamples)
	allGroups <- unique( myAnnT[[ groupColumn]])
	Ngrps <- length( allGroups)

	setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	# make the paths to where we read gene lists and write results
	optT <- readOptionsTable( optionsFile)
	if (is.null( results.path)) results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	dePath <- file.path( results.path, toolName, paste( prefix, folderName, sep="."))
	outPath <- file.path( dePath, descriptor)
	if ( ! file.exists( outPath)) dir.create( outPath, recursive=T)
	suffix <- c( "RR", "RP","SAM","EdgeR","Meta")[ match( toolName, c("RoundRobin","RankProduct","SAM","EdgeR","MetaResults"))]

	for ( i in 1:Ngrps) {

		thisGroup <- allGroups[i]
		# get the DE file of gene ratios, that has the geneIDs in ranked order
		filein <- file.path( dePath, paste( thisGroup, prefix, suffix, "Ratio.txt", sep="."))
		tbl <- read.delim( filein, as.is=T)
		if ( ! (geneColumn %in% colnames(tbl))) {
			cat( "\nGene column not found.  Looked for: ", geneColumn, "\nFound: ", colnames(tbl))
			stop()
		}
		geneList <- tbl[[ geneColumn]]
		Ngenes <- length( geneList)

		# do about 10 iterations, covering the first 1/3 of the genome.
		stopN <- round( (Ngenes*0.33) / 100) * 100
		startN <- round( (stopN*0.1) / 50) * 50
		stepN <- startN
		steps <- seq( startN, stopN, by=stepN)

		# do those enrichment calls
		ans <- geneSetMetaEnrichment( geneList, Ngenes=steps, geneSets=geneSets, wt.enrich=wt.enrich, 
					wt.pvalue=wt.pvalue, verbose=F)
		ans <- subset( ans, AVG_PVALUE <= maxPvalue)
		rownames(ans) <- 1:nrow(ans)

		# turn it to a HTML result
		outfile <- file.path( outPath, paste( thisGroup, "Enrichment.html", sep="."))
		title <- paste( "Pathway Enrichment: &nbsp; Meta Results: &nbsp; Up in '", thisGroup, "'", 
				" &nbsp;  Species: ", speciesID, sep="")
		metaEnrichment2html( ans, file=outfile, title=title)
		cat( "\nWrote file: ", outfile)
	}
	cat( "\nDone.\n")

	if (verbose) {
		cat( "\n\n=============================")
		cat( "\n\nFinished pipe 'GeneSetEnrichment' on Sample Set:     ", unlist(sampleIDset), 
			"\nSpecies: ", speciesID,"\n")
	}

	return()
}
