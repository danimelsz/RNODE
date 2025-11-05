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
| *sharedNodes*             | Compare support values of shared clades between two trees. The outputs are (1) basic statistics about number of shared clades and support values; (2) a dataframe with node labels, descendants, and support values of shared clades, which facilitates descriptive and statistical comparisons of clade composition and support between corresponding nodes.  |
| *uniqueNodes*             | Data 3    |
| *retrodictNodes*          | Data 3    |
| *normalizedSPR*           | Data 3    |
| *multiSPR*                | Data 3    |
| *summaryTopologicalDist*  | Data 3    |
| *filterMissing*           | Data 3    |
| *splitOrdFromUnord*       | Data 3    |
| *mapSupport*              | Data 3    |

The following examples are designed for users with little experience. If you have questions, send a message using GitHub issues. 

### Example 1: Comparison of support values between trees

### Example 2: Comparison 

## Cite

If you use **RNODE**, please cite this repository.
