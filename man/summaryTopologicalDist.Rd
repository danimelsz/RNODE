% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summaryTopologicalDist.R
\name{summaryTopologicalDist}
\alias{summaryTopologicalDist}
\title{summaryTopologicalDist}
\usage{
summaryTopologicalDist(
  tree1,
  tree2,
  shared = T,
  unique = T,
  RF = T,
  CID,
  outgroup = NULL,
  root = NULL
)
}
\arguments{
\item{tree1}{A .phylo tree that can be loaded using ape::read.tree for NEWICK files or TreeTools::ReadTntTree for TNT files}

\item{tree2}{Another .phylo tree}

\item{shared}{Optional. Number of shared clades between trees (default: T)}

\item{unique}{Optional. Number of unique clades in tree 1 and 2 (default: T)}

\item{RF}{Optional. Calculates the RF distance between unrooted trees based on Penny and Hendy (1985) (default: T)}

\item{CID}{Optional. Calculates the CID between trees based on Smith (2020) (default: T)}
}
\description{
suummaryTopDistances() summarizes metrics of topological distances.
}
\examples{
# Example 1 (Calculates all metrics)
tree1 = read.tree (text="(t1,(t2,(t3,(t4,t5)75)32)45);")
tree2 = read.tree (text="(t1,(t6,(t3,(t4,t5)47)53)94);")
summaryTopDistances (tree1, tree2)

# Example 2 (Calculates all metrics except the number of shared clades)
tree1 = read.tree (text="(t1,(t2,(t3,(t4,t5)75)32)45);")
tree2 = read.tree (text="(t1,(t6,(t3,(t4,t5)47)53)94);")
summaryTopDistances (tree1, tree2, shared = F)

}
\author{
Daniel YM Nakamura
}
