\name{geneSetEnrichment}
\alias{geneSetEnrichment}
\title{
Gene Pathway Enrichment
}
\description{
Find the pathways enriched for in a set of genes.
}
\usage{
geneSetEnrichment(genes, geneSets = c("GO.BiologicalProcess", "GO.MolecularFunction", 
		"GO.CellularComponent", "KEGG.Pathways", "MetabolicPathways", "GeneProduct"), 
		upOnly = TRUE, minEnrich = 2, maxPvalue = 0.05, geneUniverse = NULL)
}
\arguments{
  \item{genes}{
a vector of genes (GeneIDs for most parasite species,  Gene Symbols for most mammals),
that represent the genes of interest, perhaps from differential expression analysis.
}
  \item{geneSets}{
the named pathway definitions as a character vector, or a list of vectors of genes.
}
  \item{upOnly}{
boolean, report only over-enriched pathways or both up and down.
}
  \item{minEnrich}{
minimum enrichment to report.  Higher value reports less pathways.
}
  \item{maxPvalue}{
maximum P-value to report.  Lower value reports less pathways.
}
  \item{geneUniverse}{
the universe of potentially detected genes.  Defaults to all the genes in the geneMap for
the current species.  Useful when only a subset of genes are eligble.
}
}
\details{
Based on the hypergeometric function.  See \code{\link{enrichment}}.  Evaluates every 
named set of genes in \code{geneSets} to see how many of the elements of \code{genes} are 
present, and if that number is greater or less than expected by chance.  Tabulated for all
and then trimmed to the subset meeting the reporting criteria.
}
\value{
Two data frames, one called \code{pathways} and one called \code{genes}.  
The pathway object is the enriched pathway groups, sorted by P-value.  The most
unexpected pathway groups will be at the top.  The genes object is all the genes 
present in any enriched pathway, with all its pathways memberships, sorted by GeneID.
The pathway data frame:
  \item{PathName }{name of the pathway}
  \item{N_Total }{total number of genes (given the geneUniverse) in that pathway}
  \item{Pct_Total }{percentage of all genes that are in that pathway}
  \item{N_Given }{number of given genes that were found in that pathway}
  \item{Pct_Given }{percentage of given genes that were found in that pathway}
  \item{Enrichment }{the enrichment is the ratio of \code{Pct_Given/Pct_Total}}
  \item{P_value }{the likelihood of getting that many given genes in a pathway of that
  size, by chance}
}
\author{
  Bob Morrison
}
\note{
  The default pathway objects are distributed as part of the DuffyTools data directory.
  Not all defined species have all types of pathway files at this time.
}

\seealso{
 \code{\link{pipe.GeneSetEnrichment}}, for a higher level wrapper for evaluating 
 enrichment of multiple differential expression results.
 \code{\link{pipe.GeneSetAnalysis}}, for differential expression of pathways.
}
\examples{
  setCurrentSpecies( "Hs_grc")
  geneSetEnrichment( c("ACVR1","BAMBI","CDK8","COL1A2","ENG","FSHB","GDF9","JUNB","LTBP1","MYC","SKI","SMAD5","TGFB1","UBB"), maxP=0.0001)
}
