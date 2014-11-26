\name{expressionFileSetToMatrix}
\alias{expressionFileSetToMatrix}
\title{
Turn Expression Files into a Matrix of Gene Values
}
\description{
Combine a set of transcriptome text files into one numeric matrix, with 
a row for each gene and a column for each dataset.
}
\usage{
expressionFileSetToMatrix(fnames, fids, geneColumn = c("GENE_ID","GeneID"),
                          intensityColumn = c("INTENSITY", "RPKM_M", "RANK"),
                          missingGenes = c("drop","fill","na"), sep="\t",
			  keepIntergenics = FALSE, verbose = FALSE)
}
\arguments{
  \item{fnames}{
character vector of full pathnames to the transcript files
}
  \item{fids}{
character vector of SampleID terms (same length as \code{fnames}) to 
become the column names of the matrix
}
  \item{geneColumn}{
the column name of the gene identifier column in all transcript files
}
  \item{intensityColumn}{
the column name of the expression magnitude column in all transcript files
}
  \item{missingGenes}{
how to handle genes that are not present in every file
}
  \item{sep}{
passed to \code{read.delim} for reading in the expression files
}
  \item{keepIntergenics}{
logical, whether to keep or drop the non-gene rows from all transcriptome files
}
}
\details{
Genes can be in any row order in the input files, and will be output in alphabetical order.  
An error occurs if any file does not exist, and the list of missing files is reported.
The handling of genes not present in every file is controllable:  'drop' removes entire rows
whenever a gene is missing from any file;  'fill' fills in missing values with the smallest
magnitude value observed;  'na' sets missing values to \code{NA}.
}
\value{
A matrix of gene expression values, having a column for each transcript file 
(with \code{fids} as column names), and a row for each gene (with GeneIDs as 
rownames).
}
\author{
Bob Morrison
}
