#' @title summaryTopologicalDist
#' @name summaryTopologicalDist
#' @description suummaryTopDistances() summarizes metrics of topological distances.
#' @author Daniel YM Nakamura
#'
#' @param tree1 A .phylo tree that can be loaded using ape::read.tree for NEWICK files or TreeTools::ReadTntTree for TNT files
#' @param tree2 Another .phylo tree
#' @param shared Optional. Number of shared clades between trees (default: T)
#' @param unique Optional. Number of unique clades in tree 1 and 2 (default: T)
#' @param RF Optional. Calculates the RF distance between unrooted trees based on Penny and Hendy (1985) (default: T)
#' @param CID Optional. Calculates the CID between trees based on Smith (2020) (default: T)
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
                               shared=T, unique=T, RF=T, CID,
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

  # Number of shared clades
  sharedClades = suppressWarnings(suppressMessages(RNODE::sharedNodes(tree1_pruned, tree2_pruned, composition=F, dataframe=F, messages=F, spearman=F)))
  nSharedClades = nrow(sharedClades)

  # Number of unique clades
  uniqueClades = uniqueNodes(tree1_pruned, tree2_pruned, composition=F, dataframe=F)
  nUniqueClades1 = nrow(uniqueClades[[1]]) # unique in tree 1
  nUniqueClades2 = nrow(uniqueClades[[2]]) # unique in tree 2

  # RF distance
  unrooted_1 = unroot(tree1_pruned)
  unrooted_2 = unroot(tree2_pruned)
  RFdist = RF.dist(unrooted_1, unrooted_2, normalize = T)

  # CID
  CID = ClusteringInfoDistance(tree1_pruned, tree2_pruned, normalize = F)

  # SPR moves


  # Print
  print(paste("No. shared clades =", nSharedClades))
  print(paste("No. unique clades in tree 1 =", nUniqueClades1))
  print(paste("No. unique clades in tree 2 =", nUniqueClades2))
  print(paste("Normalized Robinson-Foulds =", RFdist))
  print(paste("Normalized Clustering Information Distance =", CID))
}

