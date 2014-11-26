\name{ortholog}
\alias{ortholog}
\title{
Gene Orthologs
}
\description{
Convert GeneIDs between organisms
}
\usage{
ortholog(genes, from = "PF3D7", to = "PY17X")
}
\arguments{
  \item{genes}{
Character vector of GeneIDs, from the species specified by \code{from}
}
  \item{from}{
Name of the current species
}
  \item{to}{
Name of the alternate species, to convert the GeneIDs to
}
}
\details{
Currently, just plasmodium species are provided.  Orthologs are currated from PlasmoDB,
OrthoMCL, and other sources.
}
\value{
A character vector of the same length as \code{genes}, holding the GeneID in the other species.
Genes with no ortholog, and genes not found in the \code{from} species, are returned as empty strings.
}
\examples{
  ortholog( "PF3D7_0102800", from="PF3D7", to="PY17X")
}
