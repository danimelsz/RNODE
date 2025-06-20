% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summaryTopologicalDist.R
\name{summaryTopologicalDist}
\alias{summaryTopologicalDist}
\title{summaryTopologicalDist}
\usage{
summaryTopologicalDist(tree1, tree2, outgroup = NULL, root = NULL)
}
\arguments{
\item{tree1}{A .phylo tree that can be loaded using ape::read.tree for NEWICK files or TreeTools::ReadTntTree for TNT files}

\item{tree2}{Another .phylo tree}

\item{outgroup}{Optional. Specify outgroup taxa to remove (by default, the function assumes that the user does not want to remove outgroup taxa)}

\item{root}{Optional. Specify the same root for both trees, which is recommended to facilitate tree comparisons (by default, the function assumes that trees share the same root)}
}
\description{
summaryTopologicalDist() summarizes metrics of topological distances (number of shared and unique clades, normalized Robinson-Foulds, normalized CID, and mean SPR moves).
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
