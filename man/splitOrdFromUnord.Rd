% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splitOrdFromUnord.R
\name{splitOrdFromUnord}
\alias{splitOrdFromUnord}
\title{splitOrdFromUnord}
\usage{
splitOrdFromUnord(input, output_index, list_ordered, invariant = F)
}
\arguments{
\item{input}{Input file (concatenated morphological matrix in Nexus or TNT format).}

\item{output_index}{Output index (e.g. if the user specify it as "Desktop/Index", the output files will be "Desktop/Index_ORDERED.nexus" and "Desktop/Index_UNORDERED.nexus")}

\item{list_ordered}{List of ordered characters e.g. c(1, 3, 9, 13). Character numbering starts with 1 (even if input data is .TNT).}

\item{invariant}{Logical value specifying if invariant characters should be removed. Following IQTREE, invariants are columns with (1) a unique state (i.e. all cells with the same value) or (2) a unique state or ? or - (i.e. some cells are the same number, other cells are ?, and/or other cells are -)}
}
\description{
Splits a morphological matrix containing ordered and unordered characters into two matrices (ordered and unordered). This is useful to run phylogenetic analyses with MK and ORDERED models in IQTREE, especially if a concatenated matrix containing invariant characters is given as input.
}
\examples{
# Example
splitOrdFromUnord

}
\author{
Daniel YM Nakamura
}
