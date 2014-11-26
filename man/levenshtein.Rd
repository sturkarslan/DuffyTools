\name{levenshtein}
\alias{levenshtein}
\alias{levenshteinDistanceMatrix}
\title{
Levenshtein String Edit Distance
}
\description{
Measure the Levenshtein edit distance between strings
}
\usage{
levenshtein(string1, string2, case = TRUE, damerau = FALSE, map = NULL, 
	    motif = FALSE)

levenshteinDistanceMatrix(stringVector, case = TRUE, damerau = FALSE, 
	   map = NULL, motif = FALSE)
}
\arguments{
  \item{string1, string2}{
character strings to measure edit distance between
}
  \item{stringVector}{
character vector to measure all pairwise edit distances of
}
  \item{case}{
be case sensitive
}
  \item{damerau}{
Demerau enhancement, swap of adjacent letters costs 1 instead of 2
}
  \item{map}{
character vector of regular expression modifiers
}
  \item{motif}{
when true, a single letter replacement costs 2, as a 'remove then insert'
}
}
\value{
For \code{levenshtein}, a numeric scalar that is the number of character edits
needed to turn \code{string1} into \code{string2}.

For \code{levenshteinDistanceMatrix}, a square symmetric numeric matrix of size 
\code{ N = length(stringVector)}, where each element is the number of character edits
between the corresponding 2 strings of \code{stringVector}.
}
\references{
Written by:   Hans-Joerg Bibiko      bibiko@eva.mpg.de
}
