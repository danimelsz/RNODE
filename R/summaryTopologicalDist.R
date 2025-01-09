#' @title summaryTopologicalDist
#' @name summaryTopologicalDist
#' @description summaryTopologicalDist() summarizes metrics of topological distances (number of shared and unique clades, normalized Robinson-Foulds, normalized CID, and mean SPR moves).
#' @author Daniel YM Nakamura
#'
#' @param tree1 A .phylo tree that can be loaded using ape::read.tree for NEWICK files or TreeTools::ReadTntTree for TNT files
#' @param tree2 Another .phylo tree
#' @param outgroup Optional. Specify outgroup taxa to remove (by default, the function assumes that the user does not want to remove outgroup taxa)
#' @param root Optional. Specify the same root for both trees, which is recommended to facilitate tree comparisons (by default, the function assumes that trees share the same root)
#'
#' @examples
#' # Example 1 (Calculates all metrics)
#' tree1 = read.tree (text="(t1,(t2,(t3,(t4,t5)75)32)45);")
#' tree2 = read.tree (text="(t1,(t6,(t3,(t4,t5)47)53)94);")
#' summaryTopDistances (tree1, tree2)
#'
#' # Example 2 (Calculates all metrics except the number of shared clades)
#' tree1 = read.tree (text="(t1,(t2,(t3,(t4,t5)75)32)45);")
#' tree2 = read.tree (text="(t1,(t6,(t3,(t4,t5)47)53)94);")
#' summaryTopDistances (tree1, tree2, shared = F)
#'
#' @export
summaryTopologicalDist = function(tree1, tree2,
                               outgroup=NULL, root=NULL){
  # Initial warnings
  missing_params <- c()
  if (is.null(tree1)) missing_params <- c(missing_params, "tree1")
  if (is.null(tree2)) missing_params <- c(missing_params, "tree2")
  if (length(missing_params) > 0) {   # Check if there are any missing parameters and print a message
    message("The following parameters are missing: ", paste(missing_params, collapse = ", "))
  } else {
    message("All required parameters provided.") # Proceed with the main functionality if all parameters are provided
  }

  #################
  # PREPROCESSING #
  #################

  # If specified, prune outgroup terminals
  if (!is.null(outgroup)) {
    tree1 <- drop.tip(tree1, outgroup)
    tree2 <- drop.tip(tree2, outgroup)
  }

  # Prune terminals not shared by both trees
  shared_terminals <- intersect(tree1$tip.label, tree2$tip.label)
  tree1_pruned <- drop.tip(tree1, tree1$tip.label[!(tree1$tip.label %in% shared_terminals)])
  tree2_pruned <- drop.tip(tree2, tree2$tip.label[!(tree2$tip.label %in% shared_terminals)])

  # If specified, reroot both trees using the same terminal
  if (!is.null(root)) {
    tree1_pruned <- root(tree1_pruned, outgroup = root)
    tree2_pruned <- root(tree2_pruned, outgroup = root)
  }

  #########################
  # TOPOLOGICAL DISTANCES #
  #########################

  # Number of polytomies
  nPolytomies1 = howManyPolytomies(tree1_pruned)
  nPolytomies2 = howManyPolytomies(tree2_pruned)

  # Number of shared clades
  sharedClades = suppressWarnings(suppressMessages(RNODE::sharedNodes(tree1_pruned, tree2_pruned, composition=F, dataframe=F, messages=F, spearman=F)))
  nSharedClades = nrow(sharedClades)

  # Number of unique clades
  nUniqueClades1 = tree1_pruned$Nnode - nSharedClades
  nUniqueClades2 = tree2_pruned$Nnode - nSharedClades

  # RF distance
  unrooted_1 = unroot(tree1_pruned)
  unrooted_2 = unroot(tree2_pruned)
  RFdist = RF.dist(unrooted_1, unrooted_2, normalize = T)

  # CID
  CID = ClusteringInfoDistance(tree1_pruned, tree2_pruned, normalize = T)

  # Print
  print(paste("No. shared clades =", nSharedClades))
  print(paste("No. unique clades in tree 1 =", nUniqueClades1))
  print(paste("No. unique clades in tree 2 =", nUniqueClades2))
  print(paste("Noo. polytomies in tree 1 =", nPolytomies1))
  print(paste("Noo. polytomies in tree 2 =", nPolytomies2))
  print(paste("Normalized Robinson-Foulds =", RFdist))
  print(paste("Normalized Clustering Information Distance =", CID))

  return(list(nSharedClades=nSharedClades,
              nUniqueClades1=nUniqueClades1,
              nUniqueClades2=nUniqueClades2,
              RFdist=RFdist,
              CID=CID,
              nPolytomies1=nPolytomies1,
              nPolytomies2=nPolytomies2))
}

