# RNODE: Comparisons of topologies and support between phylogenetic trees

[![language](https://img.shields.io/badge/language-R-blue?style=flat&logo=r&logoColor=white)](https://www.r-project.org)
[![author](https://img.shields.io/badge/author-DYM_Nakamura-blue?logo=googlescholar&logoColor=white)](https://scholar.google.com/citations?user=c0W8Cm8AAAAJ&hl=en)
[![license](https://img.shields.io/badge/license-GPL_v3-blue?logo=gnu&logoColor=white)](https://www.gnu.org/licenses/gpl-3.0.html)

**RNODE** is an R package to facilitate pre- and postprocessing of phylogenetic analyses, including (1) comparisons of topologies, branch lengths, support values, (2) comparison of DNA sequences, (3) manipulation of cladistic matrices, and (4) manipulation of trees.

Copyright (C) Daniel Y. M. Nakamura 2025

## Installation

**RNODE** can be installed with the following command:

```
devtools::install_github("danimelsz/RNODE")
```

Dependencies are expected to be automatically installed. 

## Usage

The following functions are available in **RNODE**:

| Function                  | Description |
|:--------------------------|:------------|
| *sharedNodes*             | Given two input trees, compare shared clades. The output is (1) basic statistics about number of shared clades and support values; (2) a dataframe with node labels, descendants, and support values of shared clades, which facilitates descriptive and statistical comparisons of clade composition and support between corresponding nodes.  |
| *uniqueNodes*             | Given two input trees, identify unique clades. The output is two lists containing unique clades and support values in each tree.  |
| *retrodictNodes*          | Given two input trees, create a dataframe containing support values of one tree and clade occurrence  of another tree. |
| *normalizedSPR*           | Given two binary trees, compute the normalized SPR distance, following Ding et al. (2011). |
| *multiSPR*                | Given two sets of binary trees (e.g. MPTs), compute (normalized) SPR distances between two randomly selected trees or between all pairs of trees (summarized as mean or minimum values). |
| *summaryTopologicalDist*  | Given two sets of trees, compute the number of shared clades, number of unique clades in each tree, Robinson-Foulds, and Cluster Information distance.  |
| *filterMissing*           | Given a matrix (.nex or .tnt), delete taxa and/or characters containing only missing data (?). |
| *splitOrdFromUnord*       | Given a morphological matrix (.nex or .tnt) and a list of ordered and unordered characters, split the matrix into two matrices. |
| *mapSupport*              | Given one tree with support values (e.g. majority consensus of bootstrap trees) and another tree without support values (e.g. strict consensus of optimal trees), map the support values from the former to the latter. |

The following examples are designed for users with little experience. If you have questions, send a message using GitHub issues. 

### Example 1: Comparison of support values between trees

### Example 2: Comparison 

## Cite

If you use **RNODE**, please cite this repository.
