\name{DNAtoAA}
\alias{DNAtoAA}
\alias{DNAtoAA.fast}
\alias{quickDNAtoAA}
\title{
DNA Codon to Amino Acid Peptide Conversion
}
\description{
Convert DNA sequences to amino acid peptides, in 1 or more reading frames
}
\usage{
DNAtoAA(dna, strand = "+", clipAtStop = FALSE)

DNAtoAA.fast(dna, clipAtStop = TRUE, readingFrames=1:6)

quickDNAtoAA(dnaVec)

DNAtoBestPeptide( dnaSet, clipAtStop = FALSE)
}
\arguments{
  \item{dna}{
character string of length 1, containing a sequence of nucleic acid one letter codes
}
  \item{dnaSet}{
character vector of DNA sequences
}
  \item{strand}{
coding strand that \code{dna} is from.  If '-', \code{dna} will be 
reverse complemented before conversion
}
  \item{clipAtStop}{
should conversion stop when encountering the first stop codon
}
  \item{readingFrames}{
a vector of reading frame numbers, from 1-6, where 4-6 are the 3 reading frames on the reverse strand
}
  \item{dnaVec}{
character vector of single letter nucleic acids
}
}
\details{
Facilitates converting DNA sequences to amino acid sequences.  The first form
is somewhat slow, and returns the result from all 3 reading frames.  The second
form is faster, and you choose which reading frames to use.  The third form is fastest,
expecting a vector of single nucleotides in frame.

All use the codon definitions table from the \code{\link{getCodonMap}}

Function \code{DNAtoBestPeptide} is a wrapper function that takes 1 or more DNA sequences, 
translates them in all 6 reading frames and selects the best amino acid sequence for each.
Best is defined as the reading frame yielding the minimum number of stop codons.
}
\value{
For \code{DNAtoAA()}, a character vector of length 3, giving the amino acid sequence in
each reading frame. 
For \code{DNAtoAA.fast()}, a character vector of length 1 to 6, giving the amino acid sequences. The names are the read frame numbers.
For \code{quickDNAtoAA()}, a character vector of single amino acids, of length 
\code{length(dnaVec)/3}.
The stop codon is explicitly included. Invalid DNA triplets, as from 'N's, cause '?'s in the amino acid sequence.

For \code{DNAtoBestPeptide()}, a character vector of the same length as \code{dnaSet}, with the one best
peptide sequence for each DNA sequence.
}

\seealso{
  \code{\link{getCodonMap}}, to get the conversion table of DNA triples to amino acids.
}
\examples{
  DNAtoAA( "ATGAAACGGTAGAAACCCGGG")
  DNAtoAA.fast( "ATGAAACGGTAGANACCCGGG", clipAtStop=TRUE, readingFrames=1:6)
}
