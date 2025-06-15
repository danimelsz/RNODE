#' @title normalizedSPR
#' @name normalizedSPR
#' @description \code{normalizedSPR} computes normalized SPR distances between two binary trees with the same set of leaves. If trees have polytomies, see \link{multiSPR}.
#' @author Daniel YM Nakamura, WC Wheeler, T Grant
#'
#' @param tree1 A \code{phylo} object that can be loaded using \code{ape::read.tree} for NEWICK files or \code{TreeTools::ReadTntTree} for TNT files.
#' @param tree2 Another \code{phylo} object
#' @param method Optional. Select the method for calculation of upper bound values: "traditional" or "Ding" (default). See the section Details.
#' @param outgroup Optional. Specify outgroup taxa to remove (by default, the function assumes that the user does not want to remove outgroup taxa).
#' @param root Optional. Specify the same root for both trees, which is necessary to make SPR distances meaningful (by default, the function assumes that trees share the same root).
#'
#' @details
#' The SPR distance between two trees \eqn{T_1} and \eqn{T_2} with \eqn{n} leaves is defined as the minimum number of SPR moves required to convert one tree into the other (Goloboff 2008).
#' The calculation of SPR distances is NP-hard and different heuristic procedures are available (e.g. Nakhleh et al. 2005; Beiko and Hamilton 2006; Goloboff 2008; Oliveira Martins 2008).
#'
#' SPR distance is a popular topological distance metric to measure incongruence between two trees. However, when SPR distances are compared across datasets (e.g. Torres et al. 2021), SPR distances may be coupled with the number of leaves and thus precluding statistical comparisons. Thus, normalization with values ranging between zero and one and taking into account \eqn{n} is required.
#'
#' Traditionally, the upper bound for SPR distances were \eqn{n-3}. However, Ding et al. (2011) proposed a refined upper bound for SPR distances for trees with \eqn{n>=4}:
#' \deqn{SPR_{\text{upper bound}} = n - 3 - \left(\frac{\sqrt{n - 2} - 1}{2}\right)}
#'
#' Thus, given the SPR distance (calculated with \code{TreeDist::SPRDist}) and the upper bound from Ding et al. (2011), normalized SPR distance is:
#' \deqn{SPR_{\text{normalized}} = \frac{SPR_{\text{distance}}}{SPR_{\text{upper bound}}}}
#'
#' @references Beiko, R.G., Hamilton, N., 2006. Phylogenetic identification of lateral genetic transfer events. BMC Evol Biol. 6, 15.
#' @references Ding, Y., Grünewald, S., Humphries, P.J., 2011. On agreement forests. J. Comb. Theory Ser. 118(7), 2059–2065.
#' @references Nakhleh, L., Ruths, D., Wang, L.-S., 2005. RIATA-HGT: a fast and accurate heuristic for reconstructing horizontal gene transfer. In: Proceedings of the Eleventh International Computing and Com- binatorics Conference (COCOON 05). Lecture Notes in Computer Science. LNCS no. 3595, Springer. pp. 84–93.
#' @references Goloboff, P.A., 2008. Calculating SPR distances between trees. Cladistics 24(4), 591–597.
#' @references Oliveira Martins, L., Leal, E., Kishino, H., 2008. Phylogenetic detection of recombination with a Bayesian prior on the distance between trees. PLoS One 3(7), e2651.
#' @references Torres, A., Goloboff, P.A., Catalano, S.A., 2021. Assessing topological congruence among concatenation-based phylogenomic approaches in empirical datasets. Mol. Phylogenet. Evol. 161, 107086.
#' @export
normalizedSPR = function(tree1, tree2,
                         method="Ding",
                         outgroup=NULL, root=NULL){

  #################
  # PREPROCESSING #
  #################

  # If specified, prune outgroup terminals
  if (!is.null(outgroup)) {
    tree1 <- drop.tip(tree1, outgroup)
    tree2 <- drop.tip(tree2, outgroup)
  }

  # If specified, reroot both trees using the same terminal
  if (!is.null(root)) {
    tree1 <- root(tree1, outgroup = root)
    tree2 <- root(tree2, outgroup = root)
  }

  # If a multiPhylo is provided, check if all trees are binary
  if (class(tree1)=="multiPhylo"){
    tree1 = tree1[sapply(tree1, is.binary)]
  }
  if (class(tree2)=="multiPhylo"){
    tree2 = tree2[sapply(tree2, is.binary)]
  }

  ###########################
  # NORMALIZED SPR DISTANCE #
  ###########################

  # 1. If both input files are phylo
  if (class(tree1)=="phylo" && class(tree2)=="phylo"){
  # No. leaves
  n = length(tree1$tip.label)
  # No. SPR moves
  SPR = TreeDist::SPRDist(tree1, tree2)

  if (method=="traditional"){
    # Upper bound
    upper_bound = n - 3
    # Normalized SPR
    normalizedSPR = SPR/upper_bound
  }

  else if (method=="Ding"){
  # Upper bound
  upper_bound = n - 3 - ((sqrt(n - 2) - 1) / 2)
  # Normalized SPR
  normalizedSPR = SPR/upper_bound
  }

  return (normalizedSPR)
  }

  # 2. If tree1 is multiPhylo
  else if (class(tree1)=="multiPhylo"){
    # No. leaves
    n = length(tree2$tip.label)
    # No. SPR moves
    SPR = TreeDist::SPRDist(tree1, tree2)

    if (method=="traditional"){
      # Upper bound
      upper_bound = n - 3
      # Normalized SPR
      normalizedSPR = SPR/upper_bound
    }

    else if (method=="Ding"){
      # Upper bound
      upper_bound = n - 3 - ((sqrt(n - 2) - 1) / 2)
      # Normalized SPR
      normalizedSPR = SPR/upper_bound
    }

    return (normalizedSPR)
  }

  # 3. If tree2 is multiPhylo
  else if (class(tree2)=="multiPhylo"){
    # No. leaves
    n = length(tree1$tip.label)
    # No. SPR moves
    SPR = TreeDist::SPRDist(tree1, tree2)

    if (method=="traditional"){
      # Upper bound
      upper_bound = n - 3
      # Normalized SPR
      normalizedSPR = SPR/upper_bound
    }

    else if (method=="Ding"){
      # Upper bound
      upper_bound = n - 3 - ((sqrt(n - 2) - 1) / 2)
      # Normalized SPR
      normalizedSPR = SPR/upper_bound
    }

    return (normalizedSPR)
  }
}
