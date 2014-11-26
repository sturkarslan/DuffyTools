# pipe.GeneSetAnalysis.R -- investigate sets of DE genes like pathway and GO groups

`pipe.GeneSetAnalysis` <- function( sampleIDset, speciesID="Pf3D7", annotationFile="Annotation.txt",
				optionsFile="Options.txt", results.path=NULL,  folderName="", 
				toolName=c( "MetaResults", "RoundRobin", "RankProduct", "SAM", "DESeq", "EdgeR"), 
				geneMapColumn=if(speciesID %in% c("Hs_grc","MacMu","Mmu_grc")) "NAME" else "GENE_ID", 
				groupColumn="Group", colorColumn="Color",
				geneSets=c("GO.BiologicalProcess", "GO.MolecularFunction",
				"GO.CellularComponent", "KEGG.Pathways", "MetabolicPathways", 
				"GeneProduct", "PBMC.GeneModules", "Blood.GeneModules"), 
				descriptor="GeneSets", 
				minGenesPerSet=if (speciesID %in% c("Hs_grc","MacMu","Mmu_grc")) 5 else 2, 
				mode=c("combined", "separate"), cutPvalue=0.01, verbose=T)
{

	toolName <- match.arg( toolName)
	mode <- match.arg( mode)

	if (verbose) {
		cat( "\n\n=============================")
		cat( "\nStarting pipe 'GeneSetAnalysis' on Sample Set: \n")
		print(sampleIDset)
		cat(  "\nUsing results from Species:  ", speciesID)
		cat(  "\nUsing results from DE tool:  ", toolName)
		cat(  "\nGene Sets to analyze:        ", geneSets,"\n")
	}

	annT <- readAnnotationTable( annotationFile)
	allSamples <- unique( unlist( sampleIDset))
	myAnnT <- subset( annT, SampleID %in% allSamples)
	allGroups <- unique( myAnnT[[ groupColumn]])

	DE_list <- readDEgroupsData( groupIDset=allGroups, speciesID=speciesID, 
				optionsFile=optionsFile, results.path=results.path, folderName=folderName,
				toolName=toolName)

	# allow several way of giving gene sets...
	GSanswer <- gatherGeneSets( geneSets, descriptor, mode=mode)

	nbig <- GSanswer$n
	bigSetOfDescriptors <- GSanswer$descriptor
	bigSetOfGeneSets <- GSanswer$geneSets

	# map from the group names to find the colors to use...
	where <- base::match( names( DE_list), myAnnT[[ groupColumn]])
	mycolors <- myAnnT[[ colorColumn]][ where]

	for ( i in 1:nbig) {
		geneSetAnalysis( DE_list, bigSetOfGeneSets[[i]], speciesID=speciesID, 
				colorset=mycolors, optionsFile=optionsFile, results.path=results.path, 
				folderName=folderName, toolName=toolName, geneMapColumn=geneMapColumn, 
				descriptor=bigSetOfDescriptors[i], minGenesPerSet=minGenesPerSet,
				cutPvalue=cutPvalue)
	}

	if (verbose) {
		cat( "\n\n=============================")
		cat( "\n\nFinished pipe 'GeneSetAnalysis' on Sample Set:     ", unlist(sampleIDset), 
			"\nSpecies: ", speciesID,"\n")
	}

	return()
}


`matrix.GeneSetAnalysis` <- function( x, groups=colnames(x), colors=1:ncol(x), speciesID="Pf3D7", 
				optionsFile="Options.txt", results.path=NULL,  folderName="", 
				geneSets=c("GO.BiologicalProcess", "GO.MolecularFunction",
				"GO.CellularComponent", "KEGG.Pathways", "MetabolicPathways", "GeneProduct"), 
				descriptor="GeneSets", 
				minGenesPerSet=if (speciesID %in% c("Hs_grc","MacMu","Mmu_grc")) 5 else 2, 
				mode=c("combined", "separate"), cutPvalue=0.05, verbose=T)
{

	mode <- match.arg( mode)

	if (verbose) {
		cat( "\n\n=============================")
		cat( "\nStarting 'GeneSetAnalysis' on a data matrix")
		cat(  "\nStated Species:              ", speciesID)
		cat(  "\nGene Sets to analyze:        ", geneSets,"\n")
	}

	# turn the columns of X into the DE data frames for each group
	DE_list <- matrixToDEgroupsData( x, groups=groups)

	# use the group names to get the right colors
	mygrps <- names(DE_list)
	mycolors <- colors[ match( mygrps, groups)]

	# allow several way of giving gene sets...
	GSanswer <- gatherGeneSets( geneSets, descriptor, mode=mode)

	nbig <- GSanswer$n
	bigSetOfDescriptors <- GSanswer$descriptor
	bigSetOfGeneSets <- GSanswer$geneSets

	for ( i in 1:nbig) {
		geneSetAnalysis( DE_list, bigSetOfGeneSets[[i]], speciesID=speciesID, 
				colorset=mycolors, optionsFile=optionsFile, results.path=results.path, 
				folderName=folderName, toolName=NULL, geneMapColumn="GENE_ID", 
				descriptor=bigSetOfDescriptors[i], minGenesPerSet=minGenesPerSet,
				cutPvalue=cutPvalue)
	}

	if (verbose) {
		cat( "\n\n=============================")
		cat( "\n\nFinished 'GeneSetAnalysis' on data matrix:")
	}

	return()
}



`gatherGeneSets` <- function( geneSets, descriptor="GeneSets", mode=c("combined", "separate")) {

	mode <- match.arg( mode)

	bigSetOfGeneSets <- vector( "mode"="list")
	bigSetOfDescriptors <- vector()
	nbig <- 0

	if ( is.character( geneSets)) {
		# we can be given a set of data object names, that hold an object named 'allGeneSets'
		for( nam in geneSets) {
			file <- paste( getCurrentSpeciesFilePrefix(), nam, sep=".")
			allGeneSets <- NULL
			data( list=file, envir=environment())
			if ( !is.null(allGeneSets)) {
				nbig <- nbig + 1
				bigSetOfDescriptors[nbig] <- nam
				bigSetOfGeneSets[[nbig]] <- allGeneSets
			}
		}
	} else {
		# old way is one explicit list
		nbig <- 1
		bigSetOfDescriptors[1] <- descriptor
		bigSetOfGeneSets[[1]] <- geneSets
	}

	# combine the separated geneSets if we want to...
	if ( mode == "combined" && nbig > 1) {
		newBig <- vector( mode="list")
		for ( i in 1:nbig) {
			thisName <- bigSetOfDescriptors[i]
			thisSet <- bigSetOfGeneSets[[i]]
			myNames <- names( thisSet)
			if ( nbig > 1) {
				myNewNames <- paste( thisName, myNames, sep=": ")
			} else {
				myNewNames <- myNames
			}
			names(thisSet) <- myNewNames
			newBig <- c( newBig, thisSet)
		}
		nbig <- 1
		bigSetOfDescriptors[1] <- "CombinedGeneSets"
		bigSetOfGeneSets[[1]] <- newBig
	}

	out <- list( "n"=nbig, "descriptor"=bigSetOfDescriptors, "geneSets"=bigSetOfGeneSets)
	return( out)
}
